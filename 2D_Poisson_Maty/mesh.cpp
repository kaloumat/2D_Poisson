#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <cmath>
#include "mesh.h"
#include <iostream>

using namespace std;
/* ----------------------------------------------------------------------------------- */
void Mesh_Allocate(mesh *p){
    // funkce pro alokaci pameti
    p->x = new double[p->nbNods]; // (double*)malloc(p->nbNods * sizeof(double)); // malloc vycleni v pocitacove pameti misto o velikost p->nbNods bytu
    p->y = new double[p->nbNods]; // (double *)malloc(p->nbNods * sizeof(double));

    p->BndrA = new int[p->nbBndrEdges]; // (int*)malloc(p->nbBndrEdges * sizeof(int));
    p->BndrB = new int[p->nbBndrEdges]; //(int *)malloc(p->nbBndrEdges * sizeof(int));
    p->BndrMark = new int[p->nbBndrEdges]; //(int *)malloc(p->nbBndrEdges * sizeof(int));

    p->TriA = new int[p->nbTriangles]; //(int *)malloc(p->nbTriangles * sizeof(int));
    p->TriB = new int[p->nbTriangles]; //(int *)malloc(p->nbTriangles * sizeof(int));
    p->TriC = new int[p->nbTriangles]; //(int *)malloc(p->nbTriangles * sizeof(int));
    p->TriMark = new int[p->nbTriangles]; //(int *)malloc(p->nbTriangles * sizeof(int));

    p->isDirichlet = new int[p->nbNods](); // (int*)calloc(p->nbNods, sizeof(int)); // calloc totez jako malloc, akorat vyclenou pamet jeste vynuluje (malloc tam necha nahodne hodnoty)
}
/* ----------------------------------------------------------------------------------- */
void Mesh_Free(mesh *p){
    // funkce pro dealokaci pameti
    delete(p->x);
    delete(p->y);

    delete(p->BndrA);
    delete(p->BndrB);
    delete(p->BndrMark);

    delete(p->TriA);
    delete(p->TriB);
    delete(p->TriC);
    delete(p->TriMark);
}
/* ----------------------------------------------------------------------------------- */
void Mesh_Load1(mesh *p, const char *fname){
    // inicializace pomocne struktury typu mesh
    mesh grid;

    // pomocnou strukturu se jmenem grid predame spolecne s "soubor.msh" funkci Mesh_Read
    Mesh_Read(&grid, fname);

    /*
        Z pomocne struktury predame odpovidaji cisla (pocty uzlu, usecek a trojuhelniku)
        do nasi struktury zadefinovane v main.cpp na kterou odkazuje ukazetelem p
        a pote pomoci techto cisel alokujeme pole o prislusnych velikostech pro souradnice,
        indexy uzlu hranicnich usecek, indexy uzlu trojuhelniku atd.
    */
    p->nbNods         = grid.nbNods;
    p->nbBndrEdges    = grid.nbBndrEdges;
    p->nbTriangles    = grid.nbTriangles;
    Mesh_Allocate(p);

    //vypis do terminalu pro kontrolu
    cout << "pocet uzlu = % d", p->nbNods << endl;
    cout << "\npocet usecek = %d", p->nbBndrEdges << endl ;
    cout << "\npocet trojuhelniku = %d\n", p->nbTriangles << endl;
}
/* ----------------------------------------------------------------------------------- */
void Mesh_Read(mesh *p, const char *fname){
    int      i, nbGeomElements; // pomocna neznama nbGeomElements, viz dale
    FILE     *fid;  // syntaxe pro praci se souborem, soubor pojmenuje interne fid (muze byt jakekoliv jmeno)
    char     buf[BUFFER_SIZE + 1]; //inicializace pole znaku o velikosti BUFFER_SIZE + 1, tedy 501 (plus jedna protoze pry nekdy pise chybu)
                                   // toto budeme pouzivat jako pole, do ktereho budeme vzdy ukladat vsechny znaky z daneho radku
    gmshline current; // inicializace struktury typu gmshline (definice v mesh.h)

    // vynulovani hledanych promenych
    p->nbNods         = 0;
    p->nbBndrEdges    = 0;
    p->nbTriangles    = 0;

    fid = fopen(fname, "r"); // otevreni souboru se jmenem fname -> "ctverec.msh" v modu r-reading (pouze nacitame hodnoty, nic nezapisujeme)
    if (fid == NULL){ // pro kontrolu, kdyby se nahodou z nejakeho duvodu nepodarilo otevrit soubor
        cout <<"Error 80, nepodarilo se otevrit soubor ve formatu .MSH" << endl;
        exit(1);
    }
    while (!feof(fid)){ // cyklus while, ktery probiha dokud nedojdeme do posledniho radku souboru (eof -> end of file)
        fgets(buf, BUFFER_SIZE, fid); // fgets -> nacte radek ze souboru fid o max velikosti BUFFER_SIZE -> znaky na radku ulozi do pole buf
        if (memcmp(buf, "$Nodes", 6) == 0){ //memcmp -> porovna velikost v bytech prvnich sesti znaku v poli buf se slovem $Nodes a kdyz se velikosti shoduji vrati 0
                                            // kdyz toto nastane vime, ze pro mesh format 2.2 0 8 se na dalsi radku nachazi cislo, ktere odpovida poctu uzlu (tento skript je tedy funkcni jen pro tuhle jednu verzi ulozeni z gmsh)
            fgets(buf, BUFFER_SIZE, fid); // nacte dalsi radek
            sscanf(buf, "%d", &p->nbNods); // a cislo ktere odpovida poctu uzlu site ulozi do nasi struktury pomoci ukazatele p
        }
        if (memcmp(buf, "$Elements", 9) == 0){ // stejny postup akorat gmsh format nerozlisuje usecky a trojuhelni, vsechno to radi do "elementu", takze v dalsim postupu musime rozlisit co je co
            fgets(buf, BUFFER_SIZE, fid); //nacte dalsi radek
            sscanf(buf, "%d", &nbGeomElements); // cislo na tomto radku si ulozime do pomocne nezname nbGeomElements
            for(i = 0; i < nbGeomElements; i++){ // cyklus pres vsechny elementy
                fgets(buf, BUFFER_SIZE, fid);  // vzdy nacte radek
                GmshLine_Read(&current, buf); // a funkci GmshLine predam dany radek(buf) a pomocnou strukturu typu gmshline se jmenem "current"
                                              // tato funkce nam rozlisi jestli je na danem radku uvedena usecku nebo o trojuhelnik podle nezname etyp
                switch (current.etyp){  // porovnavam hodnotu nezname etyp
                case GMSH_SEGMENT:  // v pripade, ze etyp = 1, tedy etyp = GMSH_SEGMENT jde o usecku
                    p->nbBndrEdges++;  // a pocet usecek se o jednu navysi
                    break;
                case GMSH_TRIANGLE: // v pripade, ze etyp = 2, tedy etyp = GMSH_TRIANGLE jde o trojuhelnik
                    p->nbTriangles++; // a pocet trojuhelniku se o jeden navysi
                    break;
                default: //chybova hlaska kdyby etyp se rovnal necemu jinymu nez 1 nebo 2
                    cout <<"\n\nError[MSH]:\t Unknown element %d, type %d!\n", i + 1, current.etyp<< endl;
                    exit(1);
                    break;
                }
            }
        }
    }
    fclose(fid); // uzavreni souboru fid
}
/* ----------------------------------------------------------------------------------- */
int GmshLine_Read(gmshline *p, const char *buf){
    // pomocne nezname
    int idx, etp, tgs, tgs1; // index //typ elementu //pocet oznaceni v gmsh
    int val;
    int tags[2], pom[3];
    int nred;

    // z daneho radku priradim prvni tri hodnoty do neznamych idx, etp, tgs
    nred = sscanf(buf, "%d %d %d", &idx, &etp, &tgs);
        // nejvice je pro nas dulezita neznama etp, ktera se v tom gmsh formatu nachazi na druhem miste a ma hodnotu
            // 1 kdyz jde o usecku
            // 2 kdyz jde o trojuhelnik

    if (nred != 3){ // v pripade ze na prvnich trech mistech nejsou ciselne hodnoty -> chybova hlaska
        cout <<"Error[MSH]:\t GMSH format is wrong, wrong line or buffer!\n"<< endl;
        exit(1);
    }

    p->etyp = etp; // do pomocne struktury gmshline ukladame etyp, abysme ho mohli pouzivat ve funkci Mesh_Read

    tags[0] = tags[1] = pom[0] = pom[1] = pom[2] = 0; // do pom[0], pom[1], pom[2] ukladame indexy uzlu usecek/trojuhelniku

    val  = tgs + etp * 10;
    nred = 0;
    switch (val){ // porovnavam hodnotu nezname val
    case 11: // v pripade ze val = 11 jde o primku -> dva vrcholy
        nred = sscanf(buf, "%d %d %d %d %d %d"   , &idx, &etp, &tgs1, &tags[0], &pom[0], &pom[1]);
        break;
    case 12: // v pripade ze val = 12 jde o primku -> dva vrcholy
        nred = sscanf(buf, "%d %d %d %d %d %d %d", &idx, &etp, &tgs1, &tags[0], &tags[1], &pom[0], &pom[1]);
        break;
    case 21: // v pripade ze val = 21 jde o trojuhelnik -> tri vrcholy
        nred = sscanf(buf, "%d %d %d %d %d %d %d", &idx, &etp, &tgs1, &tags[0], &pom[0], &pom[1], &pom[2]);
        break;
    case 22: // v pripade ze val = 22 jde o trojuhelnik -> tri vrcholy
        nred = sscanf(buf, "%d %d %d %d %d %d %d %d", &idx, &etp, &tgs1, &tags[0], &tags[1], &pom[0], &pom[1], &pom[2]);
        break;
    default: // chybova hlaska
        cout << "Error[MSH]:\t Incorrect number of tags (NTAGS %d)\n", tgs<< endl;
        break;
    }

    p->markPhysical    = tags[0]; // ukladam specificke cislo pro usecku/trojuhelnik
    switch (etp){ // porovnavam hodnotu nezname etp
    case GMSH_SEGMENT: // v pripade, ze etp = 1, tedy etp = GMSH_SEGMENT jde o usecku
        if (nred != tgs + 3 + 2) // chybova hlaska jestli pocet ciselnych hodnot neodpovida
            cout <<"Error[MSH]:\t Incorrect mesh format.\n"<< endl;
        p->ilist[0] = pom[0] - 1; // v pomocne strukture typu gmshline vyplnim pomoci ukazatele p pole ilist[0]
                                        // ilist[0] odpovida indexu prvniho uzlu, ktery je ulozen v pom[0] (index je snizen o jedna kvuli ccku)

        p->ilist[1] = pom[1] - 1; // v pomocne strukture typu gmshline vyplnim pomoci ukazatele p pole ilist[1]
                                        // ilist[1] odpovida indexu druheho uzlu, ktery je ulozen v pom[1] (index je snizen o jedna kvuli ccku)
        break;
    case GMSH_TRIANGLE: // v pripade, ze etp = 2, tedy etp = GMSH_TRIANGLE jde o trojuhelnik
        if (nred != tgs + 3 + 3) // chybova hlaska jestli pocet ciselnych hodnot neodpovida
            cout <<"Error[MSH]:\t Incorrect mesh format.\n"<< endl;
        p->ilist[0] = pom[0] - 1; // v pomocne strukture typu gmshline vyplnim pomoci ukazatele p pole ilist[0]
                                        // ilist[0] odpovida indexu prvniho uzlu, ktery je ulozen v pom[0] (index je snizen o jedna kvuli ccku)
        p->ilist[1] = pom[1] - 1; // v pomocne strukture typu gmshline vyplnim pomoci ukazatele p pole ilist[1]
                                        // ilist[1] odpovida indexu druheho uzlu, ktery je ulozen v pom[1] (index je snizen o jedna kvuli ccku)
        p->ilist[2] = pom[2] - 1; // v pomocne strukture typu gmshline vyplnim pomoci ukazatele p pole ilist[2]
                                        // ilist[2] odpovida indexu tretiho uzlu, ktery je ulozen v pom[2] (index je snizen o jedna kvuli ccku)
        break;
    }
    return 0;

}
/* ----------------------------------------------------------------------------------- */
void Mesh_Load2(mesh *p, const char *fname){
    // tato funkce uz jen prirazuje konkretni hodnoty do jiz alokovanych poli
    int      i, nbGeomElements, nepotrubuju[p->nbNods]; //pomocne nezname
    int      i1, i2; i1 = i2 = 0;
    FILE     *fid;
    char     buf[BUFFER_SIZE + 1];
    gmshline current;

    fid = fopen(fname, "r"); // otevru znova soubor.msh v modu r - reading
    if (fid == NULL){ //chybova hlaska
        cout <<"Error 193, nepodarilo se otevrit soubor ve formatu .MSH"<< endl;
        exit(1);
    }
    while (!feof(fid)){ // soubor ctu az do konce eof - end of file
        fgets(buf, BUFFER_SIZE, fid); // opet hledam cast kde se nachazi "$Nodes"
        if (memcmp(buf, "$Nodes", 6) == 0){
            fgets(buf, BUFFER_SIZE, fid); // na dalsim radku se nachazi pocet uzlu, to jiz nepotrebuji
            for (i = 0; i < p->nbNods; i++){
                fgets(buf, BUFFER_SIZE, fid);
                sscanf(buf, "%d %lf %lf", &nepotrubuju[i], &p->x[i], &p->y[i]); // jelikoz mam jiz alokovanou pamet mohu do pole x a y ve strukture typu mesh jiz ukladat souradnice x,y jednotlivych uzlu site
            }
        }
        if (memcmp(buf, "$Elements", 9) == 0){ // opet hledam radek, kde je napsana "$Elements"
            fgets(buf, BUFFER_SIZE, fid); // dalsi radek
            sscanf(buf, "%d", &nbGeomElements); // ulozim si pomocnou promenou nbGeomElements
            for(i = 0; i < nbGeomElements; i++){
                fgets(buf, BUFFER_SIZE, fid);
                GmshLine_Read(&current, buf); // opet pouzivam funkci pro rozliseni jestli jde o radek s informace o usecce nebo o trojuhelniku
                switch (current.etyp){ // porovnavam neznamou etyp
                case GMSH_SEGMENT: // v pripade, ze etyp = 1, tedy etyp = GMSH_SEGMENT jde o usecku
                    p->BndrA[i1] = current.ilist[0]; // ukladam do pole BndrA ve strukture typu mesh index uzlu jednoho vrcholu hranicni usecky
                    p->BndrB[i1] = current.ilist[1]; // ukladam do pole BndrB ve strukture typu mesh index uzlu druheho vrcholu hranicni usecky
                    p->BndrMark[i1] = current.markPhysical; // ukladam do pole BndrMark ve strukture typu mesh specificke cislo pro usecku
                    if(current.markPhysical == 111) // jestli je na usecce predepsan dirichlet
                        p->isDirichlet[p->BndrA[i1]] = 1; // tak si zapisuju do pole isDirichlet ve strukture typu mesh k danemu indexu 1
                    i1++;
                    break;
                case GMSH_TRIANGLE: // v pripade, ze etyp = 2, tedy etyp = GMSH_SEGMENT jde o trojuhelnik
                    p->TriA[i2] = current.ilist[0]; // ukladam do pole TriA ve strukture typu mesh index uzlu vrcholu A
                    p->TriB[i2] = current.ilist[1]; // ukladam do pole TriB ve strukture typu mesh index uzlu vrcholu B
                    p->TriC[i2] = current.ilist[2]; // ukladam do pole TriC ve strukture typu mesh index uzlu vrcholu C
                    p->TriMark[i2] = current.markPhysical; // ukladam do pole TriMark ve strukture typu mesh specificke cislo pro trojuhelnik
                    i2++;
                    break;
                default: // chybova hlaska
                    cout <<"\n\nError[MSH]:\t Unknown element %d, type %d!\n", i + 1, current.etyp<< endl;
                    exit(1);
                    break;
                }
            }
        }
    }
    fclose(fid);


    /*
        Slouzi pro kontrolu (soubor .msh muzes otevrit primo v textovem editoru a podivat
            se primo jak vypada ta struktura dat z gmsh. V tom .msh souboru te v ty
            strukture zacinajici $Elements zajimaji vzdycky posledni dve cisla, kdyz jde
            o hranicni usecku, nebo posledni tri cisla, kdyz jde o trojuhelnik, to jsou
            prave indexy uzlu, ktere danou hranicni usecku nebo trojuhelnik tvori)
        ---
        Bacha, indexy uzlu je oproti souboru .msh o jedna snizeny, aby se s nim v ccku
        dale lepe pracovalo, takze indexy uzlu zacinaji zde v ccku od nuly, ale kdyz se
        podivas primo do toho souboru .msh tak tam zacinaji od jednicky.
    */
    cout << "\nsouradnice x,y:\n"<< endl;
    for (i = 0; i < p->nbNods; i++)
        cout <<"%lf %lf\n", p->x[i], p->y[i]<< endl ;
    cout <<"\nhranicni usecky:\n" << endl;
    for (i = 0; i < p->nbBndrEdges; i++)
        cout <<"%d %d %d\n", p->BndrA[i], p->BndrB[i], p->BndrMark[i]<< endl;
    cout <<"\nisDirichlet:\n" << endl;
    for(i = 0; i < p->nbNods; i++)
        cout <<"uzel [%d] = %d\n", i, p->isDirichlet[i]<< endl;
    cout <<"\ntrojuhelniky:\n" << endl;
    for (i = 0; i < p->nbTriangles; i++)
        cout <<"%d %d %d %d\n", p->TriA[i], p->TriB[i], p->TriC[i], p->TriMark[i]<< endl;

}
