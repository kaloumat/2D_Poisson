#include <stdio.h>
#include "mesh.h"

int main(){
    /*
        inicializace structury typu mesh, ktera je definovana v hlavickovem
        souboru mesh.h, a jeji pojmenovani "sit" (muze mit jakykoliv nazev)
    */
    mesh sit;


    /*
        Alokace pameti pro
        pole(array), ktere tvori souradnice (x,y) jednotlivych uzlu site

        pole, ktere tvori indexy uzlu (BndrA, BndrB) jednotlivych hranicnich usecek
        pole, ktere tvori specificke cisla (BndrMark), tyto cisla oznacuji jaka okrajova
            podminka je na dane hranicni usecce predepsana, pro Dirichleta(BndrMark = 111) a
            pro Neumannova (BndrMark = 222) -> v nasem pripade je to vzdy 111
            ---
            specificke cislo muze byt jakekoliv, zalezi na nasi domluve, nastuvuje se v
            programu gmsh, pri vytvareni site

        pole, ktere tvori indexy uzlu (TriA, TriB, TriC) jednotlivych trojuhelniku site
        pole, ktere tvori specificke cisla (TriMark), ktera oznacuji v jake oblasti se
            dany trojuhelnik nachazi -> v nasem pripade je oblast pouze jedna, takze to
            nehraje roli, ale v BP jsem mel ty oblasti dve a pro kazdou jsem v pripade
            tepla nastavoval jinej soucinitel tepelny vodivosti
            ---
            opet se specificke cislo nastavuje v programu gmsh, pri vytvareni site, ja
            jsem napriklad zvolil cislo 4444

        pole, ktere tvori jednicky, v pripade, ze je v uzlu predepsana Dirichletova okr.
            podminka a nuly kde zadna okrajova podminka neni (mohla by jeste obsahovat
            napr dvojky, pro pripad Neumannovi okr. podminky)
            ---
            nam nestaci znat pouze hranicni usecky na kterych jsou okrajove podminky predepsany,
            ale musime to vedet primo o danych uzlech, ktere tyto usecky tvori, protoze vsechny
            vypocty probihaji prave v uzlech site
        VIZ Mesh_Allocate v mesh.cpp
        ---
        Pomoci Mesh_Load1 tak nejprve zjistime z uvedeneho souboru ("ctverect.txt")
        dulezita cisla a to
            nbNods - pocet uzlu
            nbBndrEdges - pocet hranicnich usecek
            nbTriangles - pocet trojuhelniku
        a podle jejich velikosti alokujeme uvedena pole.
    */
    Mesh_Load1(&sit, "ctverecnxn.msh");


    /*
        Kdyz uz mame alokovana pole o prislusne velikost, tak pomoci Mesh_Load2 do nich
        zapiseme ze souboru "ctverec.msh" odpovidajici hodnoty tedy
            do pole (x,y) zapiseme konkretni souradnice
            do pole (BndrA, BndrB) zapiseme konkretni indexy uzlu danych hranicnich usecek
            atd.
    */
    Mesh_Load2(&sit, "ctverecnxn.msh");


    /*
        Kdyz uz s uvedenymi poli nepotrebujeme pracovat dealokujeme pamet (v kompletnim
        skriptu bude prikaz uveden az na konci), viz funkce Mesh_Free v mesh.cpp
    */
    Mesh_Free(&sit);

    return 0;
}
