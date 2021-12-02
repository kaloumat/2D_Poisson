#define DEBUG
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <cmath>
#include "mesh.h"
#include <iostream>

using namespace std;
/* ----------------------------------------------------------------------------------- */
void Mesh::Allocate(){
    // funkce pro alokaci pameti
    x = new double[nbNods];
    y = new double[nbNods];

    BndrA = new int[nbBndrEdges];
    BndrB = new int[nbBndrEdges];
    BndrMark = new int[nbBndrEdges];

    TriA = new int[nbTriangles];
    TriB = new int[nbTriangles];
    TriC = new int[nbTriangles];
    TriMark = new int[nbTriangles];

    isDirichlet = new int[nbNods]();
}
/* ----------------------------------------------------------------------------------- */
void Mesh::Free(){
    // funkce pro dealokaci pameti
    delete(x);
    delete(y);

    delete(BndrA);
    delete(BndrB);
    delete(BndrMark);

    delete(TriA);
    delete(TriB);
    delete(TriC);
    delete(TriMark);

    delete(isDirichlet);
}
/* ----------------------------------------------------------------------------------- */
void Mesh::Load(const char *fname){
    // inicializace pomocne tridy typu Mesh
    Mesh grid;

    // pomocne classe predame "soubor.msh" a pouzijeme metodu(funkce ve tride Mesh) .Read
    grid.Read(fname);

    /*
        Z pomocne tridy "grid" predame odpovidaji atributy (pocty uzlu, usecek a trojuhelniku)
        do nasi hlavni tridy zadefinovane v main.cpp a pote pomoci techto cisel alokujeme pole
        o prislusnych velikostech pro souradnice, indexy uzlu hranicnich usecek, indexy uzlu trojuhelniku atd.
    */
    nbNods      = grid.nbNods;
    nbBndrEdges = grid.nbBndrEdges;
    nbTriangles = grid.nbTriangles;
    Allocate();

    // vypis do terminalu pro kontrolu
#ifdef DEBUG
    /*
        Kdyz v uvodu tohoto souboru zakomentujeme "#define DEBUG" tak se cast kodu ohranicena
        vyrazy "#ifdef DEBUG" a "#endif" nespusti -> nemusime jednotlive zakomentovat.
    */
    printf("pocet uzlu = %d", nbNods ); //COUT MI PISE CHYBU :/
/*    cout << "pocet uzlu = %d", this->nbNods << endl;
    cout << "\npocet usecek = %d", nbBndrEdges << endl ;
    cout << "\npocet trojuhelniku = %d\n", nbTriangles << endl; */
#endif // DEBUG

    // !!!!!DALE JSEM KOPIRUJU CAST Z KODU Z PUVODNI FUNKCE MESH_LOAD2!!!!!
    // hlavni vyhodu, ze nemusime mit Mesh_Load1 a Mesh_Load2 a staci jen Mesh_Load
    // tato funkce uz jen prirazuje konkretni hodnoty do jiz alokovanych poli

    int      i, nbGeomElements, nepotrubuju; //pomocne nezname
    int      i1, i2; i1 = i2 = 0;
    FILE     *fid;
    char     radek[500]; //prejmenoval jsem na "buf" na "radek" a "BUFFER_SIZE" jsem vymazal z globalnich velicin
    gmshline current; // inicializace pomocne tridy typu gmshline

    fid = fopen(fname, "r"); // otevru znova soubor.msh v modu r - reading
    if (fid == NULL){ //chybova hlaska
        cout <<"Error, nepodarilo se otevrit soubor ve formatu .MSH"<< endl;
        exit(1);
    }
    while (!feof(fid)){ // soubor ctu az do konce eof - end of file
        fgets(radek, sizeof(radek), fid); // opet hledam cast kde se nachazi "$Nodes"
        if (memcmp(radek, "$Nodes", 6) == 0){
            fgets(radek, sizeof(radek), fid); // na dalsim radku se nachazi pocet uzlu, to jiz nepotrebuji
            for (i = 0; i < nbNods; i++){
                fgets(radek, sizeof(radek), fid);
                sscanf(radek, "%d %lf %lf", &nepotrubuju, &x[i], &y[i]); // jelikoz mam jiz alokovanou pamet mohu do pole x a y ve strukture typu mesh jiz ukladat souradnice x,y jednotlivych uzlu site
            }
        }
        if (memcmp(radek, "$Elements", 9) == 0){ // opet hledam radek, kde je napsana "$Elements"
            fgets(radek, sizeof(radek), fid); // dalsi radek
            sscanf(radek, "%d", &nbGeomElements); // ulozim si pomocnou promenou nbGeomElements
            for(i = 0; i < nbGeomElements; i++){
                fgets(radek, sizeof(radek), fid);
                current.Read(radek); // opet pouzivam funkci pro rozliseni jestli jde o radek s informace o usecce nebo o trojuhelniku
                switch (current.etyp){ // porovnavam neznamou etyp
                case GMSH_SEGMENT: // v pripade, ze etyp = 1, tedy etyp = GMSH_SEGMENT jde o usecku
                    BndrA[i1] = current.ilist[0]; // ukladam do pole BndrA ve tride typu Mesh index uzlu jednoho vrcholu hranicni usecky
                    BndrB[i1] = current.ilist[1]; // ukladam do pole BndrB ve tride typu Mesh index uzlu druheho vrcholu hranicni usecky
                    BndrMark[i1] = current.markPhysical; // ukladam do pole BndrMark ve tride typu Mesh specificke cislo pro usecku
                    if(current.markPhysical == 111) // jestli je na usecce predepsan dirichlet
                        isDirichlet[BndrA[i1]] = 1; // tak si zapisuju do pole isDirichlet ve tride typu Mesh k danemu indexu 1
                    i1++;
                    break;
                case GMSH_TRIANGLE: // v pripade, ze etyp = 2, tedy etyp = GMSH_SEGMENT jde o trojuhelnik
                    TriA[i2] = current.ilist[0]; // ukladam do pole TriA ve tride typu Mesh index uzlu vrcholu A
                    TriB[i2] = current.ilist[1]; // ukladam do pole TriB ve tride typu Mesh index uzlu vrcholu B
                    TriC[i2] = current.ilist[2]; // ukladam do pole TriC ve tride typu Mesh index uzlu vrcholu C
                    TriMark[i2] = current.markPhysical; // ukladam do pole TriMark ve tride typu Mesh specificke cislo pro trojuhelnik
                    i2++;
                    break;
                default: // chybova hlaska
                    printf("\n\nError[MSH]:\t Neznamy element %d, type %d!\n", i + 1, current.etyp);
                    //cout <<"\n\nError[MSH]:\t Neznamy element %d, type %d!\n", i + 1, current.etyp<< endl; //COUT MI PISE CHYBU :/
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
    for (i = 0; i < nbNods; i++)
       // cout <<"%lf %lf\n", x[i], y[i]<< endl ; //COUT MI PISE CHYBU :/
        printf("%lf %lf\n", x[i], y[i]);
    cout <<"\nhranicni usecky:\n" << endl;
    for (i = 0; i < nbBndrEdges; i++)
        //cout <<"%d %d %d\n", BndrA[i], BndrB[i], BndrMark[i]<< endl; //COUT MI PISE CHYBU :/
        printf("%d %d %d\n", BndrA[i], BndrB[i], BndrMark[i]);
    cout <<"\nisDirichlet:\n" << endl;
    for(i = 0; i < nbNods; i++)
       // cout <<"uzel [%d] = %d\n", i, isDirichlet[i]<< endl; //COUT MI PISE CHYBU :/
       printf("uzel [%d] = %d\n", i, isDirichlet[i]);
    cout <<"\ntrojuhelniky:\n" << endl;
    for (i = 0; i < nbTriangles; i++)
        //cout <<"%d %d %d %d\n", TriA[i], TriB[i], TriC[i], TriMark[i]<< endl; //COUT MI PISE CHYBU :/
        printf("%d %d %d %d\n", TriA[i], TriB[i], TriC[i], TriMark[i]);
}
/* ----------------------------------------------------------------------------------- */
void Mesh::Read(const char *fname){
    int      i, nbGeomElements; // pomocna neznama nbGeomElements, viz dale
    FILE     *fid;  // syntaxe pro praci se souborem, soubor pojmenuje interne fid (muze byt jakekoliv jmeno)
    char     radek[500]; //prejmenoval jsem na "buf" na "radek" a "BUFFER_SIZE" jsem vymazal z globalnich velicin
                        // toto budeme pouzivat jako pole, do ktereho budeme vzdy ukladat vsechny znaky z daneho radku
    gmshline current; // inicializace pomocne tridy typu gmshline (definice v mesh.h)

    // vynulovani hledanych promenych
    nbNods      = 0;
    nbBndrEdges = 0;
    nbTriangles = 0;

    fid = fopen(fname, "r"); // otevreni souboru se jmenem fname -> "ctverec.msh" v modu r-reading (pouze nacitame hodnoty, nic nezapisujeme)
    if (fid == NULL){ // pro kontrolu, kdyby se nahodou z nejakeho duvodu nepodarilo otevrit soubor
        cout <<"Error, nepodarilo se otevrit soubor ve formatu .MSH" << endl;
        exit(1);
    }
    while (!feof(fid)){ // cyklus while, ktery probiha dokud nedojdeme do posledniho radku souboru (eof -> end of file)
        fgets(radek, sizeof(radek), fid); // fgets -> nacte radek ze souboru fid o max velikosti sizeof(radek) -> znaky na radku ulozi do pole radek
        if (memcmp(radek, "$Nodes", 6) == 0){ //memcmp -> porovna velikost v bytech prvnich sesti znaku v poli radek se slovem $Nodes a kdyz se velikosti shoduji vrati 0
                                            // kdyz toto nastane vime, ze pro mesh format 2.2 0 8 se na dalsi radku nachazi cislo, ktere odpovida poctu uzlu (tento skript je tedy funkcni jen pro tuhle jednu verzi ulozeni z gmsh)
            fgets(radek, sizeof(radek), fid); // nacte dalsi radek
            sscanf(radek, "%d", &nbNods); // a cislo ktere odpovida poctu uzlu site ulozi do nasi struktury pomoci ukazatele p
        }
        if (memcmp(radek, "$Elements", 9) == 0){ // stejny postup akorat gmsh format nerozlisuje usecky a trojuhelni, vsechno to radi do "elementu", takze v dalsim postupu musime rozlisit co je co
            fgets(radek, sizeof(radek), fid); //nacte dalsi radek
            sscanf(radek, "%d", &nbGeomElements); // cislo na tomto radku si ulozime do pomocne nezname nbGeomElements
            for(i = 0; i < nbGeomElements; i++){ // cyklus pres vsechny elementy
                fgets(radek, sizeof(radek), fid);  // vzdy nacte radek
                current.Read(radek); // a metode .Read ve tride typu gmshline predam dany radek
                                              // tato metoda nam rozlisi jestli je na danem radku uvedena usecku nebo o trojuhelnik podle nezname etyp
                switch (current.etyp){  // porovnavam hodnotu nezname etyp
                case GMSH_SEGMENT:  // v pripade, ze etyp = 1, tedy etyp = GMSH_SEGMENT jde o usecku
                    nbBndrEdges++;  // a pocet usecek se o jednu navysi
                    break;
                case GMSH_TRIANGLE: // v pripade, ze etyp = 2, tedy etyp = GMSH_TRIANGLE jde o trojuhelnik
                    nbTriangles++; // a pocet trojuhelniku se o jeden navysi
                    break;
                default: //chybova hlaska kdyby etyp se rovnal necemu jinymu nez 1 nebo 2
                    //cout <<"\n\nError[MSH]:\t Neznamy element %d, type %d!\n", i + 1, current.etyp<< endl; //COUT MI PISE CHYBU :/
                    printf("\n\nError[MSH]:\t Neznamy element %d, type %d!\n", i + 1, current.etyp);
                    exit(1);
                    break;
                }
            }
        }
    }
    fclose(fid); // uzavreni souboru fid
}
/* ----------------------------------------------------------------------------------- */
int gmshline::Read(const char *radek){
    // pomocne nezname
    int idx, etp, tgs, tgs1; // index //typ elementu //pocet oznaceni v gmsh
    int val;
    int tags[2], pom[3];
    int nred;

    // z daneho radku priradim prvni tri hodnoty do neznamych idx, etp, tgs
    nred = sscanf(radek, "%d %d %d", &idx, &etp, &tgs);
        // nejvice je pro nas dulezita neznama etp, ktera se v tom gmsh formatu nachazi na druhem miste a ma hodnotu
            // 1 kdyz jde o usecku
            // 2 kdyz jde o trojuhelnik

    if (nred != 3){ // v pripade ze na prvnich trech mistech nejsou ciselne hodnoty -> chybova hlaska
        cout <<"Error[MSH]:\t GMSH format je spatne, nebo spatny radek\n"<< endl;
        exit(1);
    }

    etyp = etp; // do pomocne struktury gmshline ukladame etyp, abysme ho mohli pouzivat ve funkci Mesh_Read
    this->idx = idx;
    tags[0] = tags[1] = pom[0] = pom[1] = pom[2] = 0; // do pom[0], pom[1], pom[2] ukladame indexy uzlu usecek/trojuhelniku

    val  = tgs + etp * 10;
    nred = 0;
    switch (val){ // porovnavam hodnotu nezname val
    case 11: // v pripade ze val = 11 jde o primku -> dva vrcholy
        nred = sscanf(radek, "%d %d %d %d %d %d"   , &idx, &etp, &tgs1, &tags[0], &pom[0], &pom[1]);
        break;
    case 12: // v pripade ze val = 12 jde o primku -> dva vrcholy
        nred = sscanf(radek, "%d %d %d %d %d %d %d", &idx, &etp, &tgs1, &tags[0], &tags[1], &pom[0], &pom[1]);
        break;
    case 21: // v pripade ze val = 21 jde o trojuhelnik -> tri vrcholy
        nred = sscanf(radek, "%d %d %d %d %d %d %d", &idx, &etp, &tgs1, &tags[0], &pom[0], &pom[1], &pom[2]);
        break;
    case 22: // v pripade ze val = 22 jde o trojuhelnik -> tri vrcholy
        nred = sscanf(radek, "%d %d %d %d %d %d %d %d", &idx, &etp, &tgs1, &tags[0], &tags[1], &pom[0], &pom[1], &pom[2]);
        break;
    default: // chybova hlaska
        //cout << "Error[MSH]:\t Spatny pocet tagu: %d\n", tgs<< endl; //COUT MI PISE CHYBU :/
        printf("Error[MSH]:\t Spatny pocet tagu: %d\n", tgs);
        break;
    }

    markPhysical = tags[0]; // ukladam specificke cislo pro usecku/trojuhelnik
    switch (etp){ // porovnavam hodnotu nezname etp
    case GMSH_SEGMENT: // v pripade, ze etp = 1, tedy etp = GMSH_SEGMENT jde o usecku
        if (nred != tgs + 3 + 2) // chybova hlaska jestli pocet ciselnych hodnot neodpovida
            cout <<"Error[MSH]:\t Spatny format site pro usecky\n"<< endl;
        ilist[0] = pom[0] - 1; // v pomocne tride typu gmshline vyplnim pole ilist[0]
                                        // ilist[0] odpovida indexu prvniho uzlu, ktery je ulozen v pom[0] (index je snizen o jedna kvuli ccku)

        ilist[1] = pom[1] - 1; // v pomocne tride typu gmshline vyplnim pole ilist[1]
                                        // ilist[1] odpovida indexu druheho uzlu, ktery je ulozen v pom[1] (index je snizen o jedna kvuli ccku)
        break;
    case GMSH_TRIANGLE: // v pripade, ze etp = 2, tedy etp = GMSH_TRIANGLE jde o trojuhelnik
        if (nred != tgs + 3 + 3) // chybova hlaska jestli pocet ciselnych hodnot neodpovida
            cout <<"Error[MSH]:\t Spatny format site pro trojuhelniky\n"<< endl;
        ilist[0] = pom[0] - 1; // v pomocne tride typu gmshline vyplnim pole ilist[0]
                                        // ilist[0] odpovida indexu prvniho uzlu, ktery je ulozen v pom[0] (index je snizen o jedna kvuli ccku)
        ilist[1] = pom[1] - 1; // v pomocne tride typu gmshline vyplnim pole ilist[1]
                                        // ilist[1] odpovida indexu druheho uzlu, ktery je ulozen v pom[1] (index je snizen o jedna kvuli ccku)
        ilist[2] = pom[2] - 1; // v pomocne tride typu gmshline vyplnim pole ilist[2]
                                        // ilist[2] odpovida indexu tretiho uzlu, ktery je ulozen v pom[2] (index je snizen o jedna kvuli ccku)
        break;
    }
    return 0;

}
