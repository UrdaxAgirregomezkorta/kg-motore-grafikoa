//	Program developed by Urdax Agirregomezkorta Beitia
//
//	Informatika Fakultatea
//	Euskal Herriko Unibertsitatea
//	http://www.ehu.eus/if
//
// to compile it: gcc dibujar-triangulos-y-objetos.c -lGL -lGLU -lglut -lm
//
//
//


#include <GL/glut.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "cargar-triangulo.h"

typedef struct mlist
    {
    double m[16];
    struct mlist *hptr;
    } mlist;

typedef struct triobj
    {
    hiruki *triptr;
    int num_triangles;
    mlist *mptr;
    struct triobj *hptr;
    unsigned char color[3]; 
    } triobj;

typedef struct eguzkiobj
    {
    double eguzkibektore[3];
    double intentsitate[3];
    } eguzkiobj;

typedef struct bonbilaobj
    {
    double puntua[3];
    double intentsitate[3];
    } bonbilaobj;



extern int load_ppm(char *file, unsigned char **bufferptr, int *dimxptr, int * dimyptr);

unsigned char *bufferra;
int dimx,dimy,dimentsioa;

int indexx;
hiruki *triangulosptr;
triobj *foptr;
triobj *sel_ptr;
int denak;
int lineak;
int objektuak;
int kamera;
char aldaketa;
int ald_lokala;
int hurbilegi;


int objektuIkuspegi;
triobj *kamptr;
int hegaldi;
double proiekMatrizea[16];
int perspektiba;
int backculling;
double Mesa[16];

int normala;
int argia;
int objektu;
char fitxiz[100];
float argiGlobala;


bonbilaobj *bonbila;
eguzkiobj *eguzkia;
int EguzkiaPiztuta;
int BonbilaPiztuta;
int ObjektufokuaPiztuta;
int KamerafokuaPiztuta;

float Ka[3];
float Kd[3];
float Ks[3];
int distira;//edo ns
double angeluaFokuaKamera;
double IKamFokua;
double IObjFokua;
double angeluaFokuaObjektua;

int zeinArgi;




void printMatrix(const char *name, double m[16]) {
    printf("%s:\n", name);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            printf("%f ", m[i + 4 * j]); // Orden de columnas
        }
        printf("\n");
    }
}
void normalize(double vec[3]) {
    double magnitude = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
    if (magnitude > 0.0) {
        vec[0] /= magnitude;
        vec[1] /= magnitude;
        vec[2] /= magnitude;
    }
}
// v1 ^ v2= RESULT
void crossProduct(const double v1[3], const double v2[3], double result[3]) {
    result[0] = v1[1] * v2[2] - v1[2] * v2[1];
    result[1] = v1[2] * v2[0] - v1[0] * v2[2];
    result[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

void lookAt(double targetX, double targetY, double targetZ) {

    double Zc[3] = { 
        (kamptr->mptr->m[12]-targetX),
        (kamptr->mptr->m[13]-targetY),
        (kamptr->mptr->m[14]-targetZ)
    };
    normalize(Zc);

    double Yc[3] = { kamptr->mptr->m[4],  kamptr->mptr->m[5], kamptr->mptr->m[6] }; 
    double Xc[3];
    crossProduct(Yc, Zc, Xc); 
    normalize(Xc);

    crossProduct(Zc, Xc, Yc);

    kamptr->mptr->m[0] = Xc[0]; kamptr->mptr->m[1] = Xc[1]; kamptr->mptr->m[2]  = Xc[2]; 
    kamptr->mptr->m[4] = Yc[0]; kamptr->mptr->m[5] = Yc[1]; kamptr->mptr->m[6]  = Yc[2];
    kamptr->mptr->m[8] = Zc[0]; kamptr->mptr->m[9] = Zc[1]; kamptr->mptr->m[10] = Zc[2];
    

}
void multMatrixVector4(double result[4], double matrix[16], double vector[4]) {

    result[0] = matrix[0] * vector[0] + matrix[4] * vector[1] + matrix[8] * vector[2] + matrix[12] * vector[3];
    result[1] = matrix[1] * vector[0] + matrix[5] * vector[1] + matrix[9] * vector[2] + matrix[13] * vector[3];
    result[2] = matrix[2] * vector[0] + matrix[6] * vector[1] + matrix[10] * vector[2] + matrix[14] * vector[3];
    result[3] = matrix[3] * vector[0] + matrix[7] * vector[1] + matrix[11] * vector[2] + matrix[15] * vector[3];
}
void multMatrixVector3(double result[3], double matrix[16], double vector[3]) {

    result[0] = matrix[0] * vector[0] + matrix[4] * vector[1] + matrix[8] * vector[2];
    result[1] = matrix[1] * vector[0] + matrix[5] * vector[1] + matrix[9] * vector[2];
    result[2] = matrix[2] * vector[0] + matrix[6] * vector[1] + matrix[10] * vector[2];
}
void multMatrixVector3Right(double result[3], double matrix[16], double vector[3]) {
    result[0] = matrix[0] * vector[0] + matrix[1] * vector[1] + matrix[2] * vector[2];
    result[1] = matrix[4] * vector[0] + matrix[5] * vector[1] + matrix[6] * vector[2];
    result[2] = matrix[8] * vector[0] + matrix[9] * vector[1] + matrix[10] * vector[2];
}
void matrizeaKalk(double modelview[16], double Mesa[16], double matrix[16]) {
    //  EMAITZA = A x B egiten du 
    for (int i = 0; i < 16; i++) {
        modelview[i] = 0.0;
    }

    // Realizar la multiplicación de matrices (formato de columnas)
    for (int i = 0; i < 4; i++) {         // Fila de la matriz resultante
        for (int j = 0; j < 4; j++) {     // Columna de la matriz resultante
            for (int k = 0; k < 4; k++) { // Columna de Mesa y fila de matrix
                modelview[i + 4 * j] += Mesa[i + 4 * k] * matrix[k + 4 * j];
            }
        }
    }
}



// biderketa kalkulatzeko
double dot_product(double v1[3], double v2[3]) {
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}


// normala kamera erreferentziara pasatzeko eta cos kalkulatzeko
double cos_angle_camera_normal(double modelview[16], hiruki *tript) {
    double normal_camera[3];
    double normal_camera_p[3];

    double triMv[3];
    double triMv_p[3];

        multMatrixVector3(normal_camera,modelview,tript->N);
        multMatrixVector3(normal_camera_p,proiekMatrizea,normal_camera);

   // printf("%f, %f, %f \n",normal_camera_p[0],normal_camera_p[1],normal_camera_p[2]);

    double tri[3] = {
     tript->p1.x,
     tript->p1.y, 
     tript->p1.z};

    multMatrixVector3(triMv,modelview,tri);

    multMatrixVector3(triMv_p,proiekMatrizea,triMv);

    double camera_direction[3] = {
     0,
     0, 
     1};

    // Calcular el coseno del ángulo entre la dirección de la cámara y el normal transformado
    double cos_theta = dot_product(camera_direction, normal_camera_p);

   // printf("%f",cos_theta );

    return cos_theta;
}


//munduaren erreferentzia sistema
void objektuari_aldaketa_sartu_ezk(double m[16])
{

    double resultado[16];
    int i, j, k;


     mlist *nueva_matriz = (mlist *)malloc(sizeof(mlist));

    for (i = 0; i < 16; i++) {
        nueva_matriz->m[i] = sel_ptr->mptr->m[i];
    }

    nueva_matriz->hptr = sel_ptr->mptr;
    sel_ptr->mptr = nueva_matriz;

    matrizeaKalk(resultado, m, sel_ptr->mptr->m);

    for (i = 0; i < 16; ++i) {
        sel_ptr->mptr->m[i] = resultado[i];
    }

    }



//objektuaren erreferentzia sistema
void objektuari_aldaketa_sartu_esk(double m[16])
{
   double resultado[16];
   double vector[3];
    int i, j, k;

   // oraingo matrizea gorde
    mlist *nueva_matriz = (mlist *)malloc(sizeof(mlist));
    if(argia==1){

        
        multMatrixVector3Right(vector,m,eguzkia->eguzkibektore);
        eguzkia->eguzkibektore[0]=vector[0];
        eguzkia->eguzkibektore[1]=vector[1];    
        eguzkia->eguzkibektore[2]=vector[2];    
        
    }else if(kamera == 1){
    for (i = 0; i < 16; i++) {
        nueva_matriz->m[i] = kamptr->mptr->m[i];  // kopiatu matrizea
    }

    nueva_matriz->hptr = kamptr->mptr;  // aurreko matrizearekin linkeatu
    kamptr->mptr = nueva_matriz;  // pointerra

    matrizeaKalk(resultado,kamptr->mptr->m,m);

    // Guardamos el resultado en la matriz actual
    for (i = 0; i < 16; ++i) {
        kamptr->mptr->m[i] = resultado[i];
    }

     }else{
    for (i = 0; i < 16; i++) {
        nueva_matriz->m[i] = sel_ptr->mptr->m[i];  // kopiatu matrizea
    }
    nueva_matriz->hptr = sel_ptr->mptr;  // aurreko matrizearekin linkeatu
    sel_ptr->mptr = nueva_matriz;  // pointerra

    matrizeaKalk(resultado,sel_ptr->mptr->m,m);


    // Guardamos el resultado en la matriz actual
    for (i = 0; i < 16; ++i) {
        sel_ptr->mptr->m[i] = resultado[i];
    }}

}
// TODO
// funtzio honek u eta v koordenatuei dagokien pointerra itzuli behar du.
// debe devolver el pointer correspondiente a las coordenadas u y v
unsigned char * color_textura(float u, float v)
{
    if (u < 0) u = 0;
    if (u > 1) u = 1;
    if (v < 0) v = 0;
    if (v > 1) v = 1;

    int desplazamendua;
    int x,y;
    unsigned char * lag;
    lag = bufferra;

    x = (int)(u * (dimx-1)) % dimx;
    y = (int)((1-v) * (dimy-1)) % dimy;

    desplazamendua = (x + y * dimx);
    return(lag + 3 * desplazamendua);
}


void kameraAnalisi(double angulo, char eje) {

  // radianetan
    double rad = angulo * 3.1415/ 180.0;

    // `at` puntura eramateko 
    double Mt_at[16] = {
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        -sel_ptr->mptr->m[12], -sel_ptr->mptr->m[13], -sel_ptr->mptr->m[14], 1
    };

   // printMatrix("Mt_at",Mt_at);
    //    printf("\n");

    // hasierako puntura joateko
    double Mtat[16] = {
         1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        sel_ptr->mptr->m[12], sel_ptr->mptr->m[13], sel_ptr->mptr->m[14], 1
    };
   //  printMatrix("Mtat",Mt_at);
     //   printf("\n");

        // rodriges matrizea...
        double MR[16];
        double c = cos(rad);    // cos(α)
        double s = sin(rad);    // sin(α)
        double t = 1 - c;       // 1 - cos(α)
        double x;
        double y;
        double z;

      //  printMatrix("kamaria biderkau baño lenuo",kamptr->mptr->m );

        if (eje == 'x') {
 
        x =kamptr->mptr->m[0];
        y =kamptr->mptr->m[1];
        z =kamptr->mptr->m[2];

   // printf("Eje Y: x = %f, y = %f, z = %f\n", x, y, z);


        } else if (eje == 'y') {
 
         x =kamptr->mptr->m[4];
         y =kamptr->mptr->m[5];
         z =kamptr->mptr->m[6];
     //   printf("Eje X: x = %f, y = %f, z = %f\n", x, y, z);
 

        }

        MR[0] = t * x * x + c;          MR[4] = t * x * y - s * z;      MR[8] = t * x * z + s * y;      MR[12] = 0;
        MR[1] = t * x * y + s * z;      MR[5] = t * y * y + c;          MR[9] = t * y * z - s * x;      MR[13] = 0;
        MR[2] = t * x * z - s * y;      MR[6] = t * y * z + s * x;      MR[10] = t * z * z + c;         MR[14] = 0;
        MR[3] = 0;                      MR[7] = 0;                      MR[11] = 0;                     MR[15] = 1;
      
      //printMatrix("MR", MR);


    double Mr1[16];
    matrizeaKalk(Mr1,MR,Mt_at);    
    // printMatrix("Mr", Mr1);
  

        double Mr2[16];

     matrizeaKalk(Mr2, Mtat, Mr1);           

    //        printMatrix("Mr", Mr2);

             double emaitza[16];
    
    matrizeaKalk(emaitza,Mr2,kamptr->mptr->m);           

                mlist *matrizeberri = (mlist *)malloc(sizeof(mlist));

                // Copiar los valores de Mr2 a matrizeberri
                     for (int i = 0; i < 16; ++i) {
                    matrizeberri->m[i] = emaitza[i];
                }

                matrizeberri->hptr = kamptr->mptr;  // Enlazar la nueva matriz a la estructura kamptr
                kamptr->mptr = matrizeberri;  // Actualizar el puntero de kamptr


                // Copiar los valores de Mr2 a la matriz de kamptr
                for (int i = 0; i < 16; ++i) {
                    kamptr->mptr->m[i] = emaitza[i];
                }



       lookAt(sel_ptr->mptr->m[12], sel_ptr->mptr->m[13], sel_ptr->mptr->m[14]);

}
void MesaLortu() {
double *M;
    if(objektuIkuspegi == 1){
    M= sel_ptr->mptr->m;
    

    }else{
    M = kamptr->mptr->m; 

    }
    double Xc[3] = { M[0], M[1], M[2] };
    double Yc[3] = { M[4], M[5], M[6] };
    double Zc[3] = { M[8], M[9], M[10] };
    double E[3] = { M[12], M[13], M[14] };

    Mesa[0] = Xc[0]; Mesa[4] = Xc[1]; Mesa[8] = Xc[2]; Mesa[12] = - (Xc[0] * E[0] + Xc[1] * E[1] + Xc[2] * E[2]);
    Mesa[1] = Yc[0]; Mesa[5] = Yc[1]; Mesa[9] = Yc[2]; Mesa[13] =  -(Yc[0] * E[0] + Yc[1] * E[1] + Yc[2] * E[2]);
    Mesa[2] = Zc[0]; Mesa[6] = Zc[1]; Mesa[10] = Zc[2]; Mesa[14] = -(Zc[0] * E[0] + Zc[1] * E[1] + Zc[2] * E[2]);
    Mesa[3] = 0.0;   Mesa[7] = 0.0;   Mesa[11] = 0.0;   Mesa[15] = 1.0;

}



void ortografikoaKalk() {
    double n = 0.0;
    double f = -100.0;
    double l = -1.0;
    double r = 1.0;
    double t = 1.0;
    double b = -1.0;

    double r_l = r - l;
    double t_b = t - b;
    double f_n = f - n;

    double r_plus_l = r + l;
    double t_plus_b = t + b;
    double f_plus_n = f + n;

    proiekMatrizea[0] = 2.0 / 2.0;
    proiekMatrizea[1] = 0.0;
    proiekMatrizea[2] = 0.0;
    proiekMatrizea[3] = 0.0;

    proiekMatrizea[4] = 0.0;
    proiekMatrizea[5] = 2.0 / 2.0;
    proiekMatrizea[6] = 0.0;
    proiekMatrizea[7] = 0.0;

    proiekMatrizea[8] = 0.0;
    proiekMatrizea[9] = 0.0;
    proiekMatrizea[10] = 2.0/100.0;
    proiekMatrizea[11] = 0.0;

    proiekMatrizea[12] =0;
    proiekMatrizea[13] = 0;
    proiekMatrizea[14] = 1.0;
    proiekMatrizea[15] = 1.0;


}

void perspektibaKalk() {

    double n = 0.1;
    double f = 100.0;
    double l = -0.1;
    double r = 0.1;
    double t = 0.1;
    double b = -0.1;


    double r_l = r - l;
    double t_b = t - b; 
    double f_n = f - n;
    double n_f= n - f;


     proiekMatrizea[0] = (2 * n) / r_l;
    proiekMatrizea[1] = 0.0;
    proiekMatrizea[2] = 0.0;
    proiekMatrizea[3] = 0.0;

    proiekMatrizea[4] = 0.0;
    proiekMatrizea[5] = (2 * n) / t_b;
    proiekMatrizea[6] = 0.0;
    proiekMatrizea[7] = 0.0;

    proiekMatrizea[8] = (r + l) / r_l;
    proiekMatrizea[9] = (t + b) / t_b;
    proiekMatrizea[10] = -100.1/(0.1-100);    
    proiekMatrizea[11] = -1.0;

    proiekMatrizea[12] = 0.0;
    proiekMatrizea[13] = 0.0;
    proiekMatrizea[14] = -2.0*0.1*100.0/(0.1-100); 
    proiekMatrizea[15] = 0.0;

    // proiekMatrizea[0] = (2 * n) / r_l;
    // proiekMatrizea[1] = 0.0;
    // proiekMatrizea[2] = 0.0;
    // proiekMatrizea[3] = 0.0;

    // proiekMatrizea[4] = 0.0;
    // proiekMatrizea[5] = (2 * n) / t_b;
    // proiekMatrizea[6] = 0.0;
    // proiekMatrizea[7] = 0.0;

    // proiekMatrizea[8] = (r + l) / r_l;
    // proiekMatrizea[9] = (t + b) / t_b;
    // proiekMatrizea[10] = -n+f / n_f;     //  -(n+f) / n_f; -(f + n) / n_f; 
    // proiekMatrizea[11] = -1.0;

    // proiekMatrizea[12] = 0.0;
    // proiekMatrizea[13] = 0.0;
    // proiekMatrizea[14] = -2 * n * f / f_n;                  //-(2 * f * n) / n_f;
    // proiekMatrizea[15] = 0.0;

 
}


void mpkalkulatu(){
if(perspektiba == 0){
ortografikoaKalk();
}else{
perspektibaKalk();
}

}


// TODO
// lerroa marrazten du, baina testuraren kodea egokitu behar da
// dimentsioa --> marrazgunearen pixel kopurua
// dibuja una linea pero hay que codificar la textura
// dimentsioa --> número de pixels de la ventana de dibujo
void  dibujar_linea_z(float linea,float c1x,float c1z, float c1u,float c1v,float c2x,float c2z,float c2u,float c2v,unsigned char color[3])
{
float xkoord,zkoord;
float cx, cz, cu, cv;
float u,v;
int i,puntukop;
unsigned char r,g,b;
unsigned char *colorv;

    r=color[0];
    g=color[1];
    b=color[2];

//printf("linea");

// TODO x balioak -1 eta 1 artekoak direla ziurtatu eta ondorioz z, u eta v egokitu

 // x balioak -1 eta 1 artekoa izan behar du, ezkerretik eskuinera 2 unitateko aldaketa dauka x balioak, eta
 // dimentsioa adina pixel behar dira ezkerretik eskuinera joateko, beraz, 2-ko aldaketari "dimentsioa" pixel dagozkio.
 // ondorioz c1x-tik c2x-ra doan tartean behar ditudan pixelak = (c2x-c1x)dimentsioa/2 pixel behar ditut
 // ondorioz, Z-ren u-ren eta v-ren aldaketan zenbaki hori erabili behar dut.
 // para un cambio de -1 a 1 en x hay que dibujar "dimentsioa" pixels, para un cambio que va de c1x a c2x, cuantos?
 puntukop = (c2x - c1x) * (float)dimentsioa / 2.0;
 //puntu bat
    if (c1x == c2x) {
        // En este caso, los puntos son el mismo, dibujar solo un punto
        glBegin(GL_POINTS);
        glColor3ub(r,g,b);
        glVertex3f(c1x, linea, c1z);
        glEnd();
        return;
    }
if (puntukop > 0) {
    // Calculate the step size for interpolation
    cx = (c2x - c1x) / (float)puntukop;
    cz = (c2z - c1z) / (float)puntukop;
    cu = (c2u - c1u) / (float)puntukop;
    cv = (c2v - c1v) / (float)puntukop;
} else {
    // If puntukop <= 0, we set the step size to 0, which stops interpolation
    cx = 0;
    cz = 0;
    cu = 0;
    cv = 0;
    }
glBegin( GL_POINTS );
    for (i = 0, xkoord=c1x, zkoord=c1z, u=c1u, v=c1v;
     i <puntukop /*xkoord <= c2x*/;
     i++,xkoord += cx /* 500 puntu -1 eta 1 artean */)
    {
    //   colorv=color_textura(u, v); // si esta función es correcta se ve la foto en la ventana
    glColor3ub(r,g,b);
    glVertex3f(xkoord, linea, zkoord );
    zkoord+=cz;
    u+=cu;
    v+=cv;
    }
glEnd();
}
// TODO (transform)
// supone que la cuarta coordenada es un 1
// honek objekturaren matrizea hartu eta objektua biderkatu behar du
// m16 erabili behar
void mxp(punto *pptr, double m[16], punto p)
{
    // Calculamos las coordenadas x', y', z' y w'
    double x_prime = m[0] * p.x + m[4] * p.y + m[8] * p.z + m[12] * 1;
    double y_prime = m[1] * p.x + m[5] * p.y + m[9] * p.z + m[13] * 1;
    double z_prime = m[2] * p.x + m[6] * p.y + m[10] * p.z + m[14] * 1;
    double w_prime = m[3] * p.x + m[7] * p.y + m[11] * p.z + m[15] * 1;

        pptr->x = x_prime / w_prime;
        pptr->y = y_prime / w_prime;
        pptr->z = z_prime / w_prime;

    pptr->u = p.u;
    pptr->v = p.v;
}


void eskala_aldaketa(double factor) {
    double eskala_matrizea[16] = {
        factor, 0,      0,      0,
        0,      factor, 0,      0,
        0,      0,      factor, 0,
        0,      0,      0,      1
    };

    if (ald_lokala) {
        objektuari_aldaketa_sartu_esk(eskala_matrizea);
    } else {
        objektuari_aldaketa_sartu_ezk(eskala_matrizea);
    }
}
void bektoreNormalaKalk(hiruki *tptr) {
    // Bi bektore kalkulatu p1-etik ateatzeianak
    double v1x = tptr->p2.x - tptr->p1.x;
    double v1y = tptr->p2.y - tptr->p1.y;
    double v1z = tptr->p2.z - tptr->p1.z;

    double v2x = tptr->p3.x - tptr->p1.x;
    double v2y = tptr->p3.y - tptr->p1.y;
    double v2z = tptr->p3.z - tptr->p1.z;

    tptr->N[0] = v1y * v2z - v1z * v2y; 
    tptr->N[1] = v1z * v2x - v1x * v2z; 
    tptr->N[2] = v1x * v2y - v1y * v2x;

    // Normalizar el vector 
    double length = sqrt(tptr->N[0] * tptr->N[0] + tptr->N[1] * tptr->N[1] + tptr->N[2] * tptr->N[2]);
    if (length > 0) {
        tptr->N[0] /= length;
        tptr->N[1] /= length;
        tptr->N[2] /= length;
    }


}
//triangelua marrazteko goiko,erdiko eta beheko puntuak aterako ditugu:
void ordenatu(punto* p1, punto* p2, punto* p3, punto** txikiena, punto** erdikoa, punto** altuena) {
if (p1->y == p2->y && p2->y == p3->y) {

        if (p1->x == p2->x && p2->x == p3->x) {
                if (p1->z >= p2->z && p1->z >= p3->z) {
                *altuena = p1;
                *erdikoa = (p2->z >= p3->z) ? p2 : p3;
                *txikiena = (p2->z >= p3->z) ? p3 : p2;
            } else if (p2->z >= p1->z && p2->z >= p3->z) {
                *altuena = p2;
                *erdikoa = (p1->z >= p3->z) ? p1 : p3;
                *txikiena = (p1->z >= p3->z) ? p3 : p1;
            } else {
                *altuena = p3;
                *erdikoa = (p1->z >= p2->z) ? p1 : p2;
                *txikiena = (p1->z >= p2->z) ? p2 : p1;
            }
        }else{
        if (p1->x >= p2->x && p1->x >= p3->x) {
            *altuena = p1;
            if (p2->x >= p3->x) {
                *erdikoa = p2;
                *txikiena = p3;
            } else {
                *erdikoa = p3;
                *txikiena = p2;
            }
        } else if (p2->x >= p1->x && p2->x >= p3->x) {
            *altuena = p2;
            if (p1->x >= p3->x) {
                *erdikoa = p1;
                *txikiena = p3;
            } else {
                *erdikoa = p3;
                *txikiena = p1;
            }
        } else {
            *altuena = p3;
            if (p1->x >= p2->x) {
                *erdikoa = p1;
                *txikiena = p2;
            } else {
                *erdikoa = p2;
                *txikiena = p1;
            }
        }
        }

     }else{
      if (p1->y >= p2->y && p1->y >= p3->y) {
        *altuena = p1;
        if (p2->y > p3->y) {
            *erdikoa = p2;
            *txikiena = p3;
        } else {
            *erdikoa = p3;
            *txikiena = p2;
        }
    } else if (p2->y >= p1->y && p2->y >= p3->y) {
        *altuena = p2;
        if (p1->y > p3->y) {
            *erdikoa = p1;
            *txikiena = p3;
        } else {
            *erdikoa = p3;
            *txikiena = p1;
        }
    } else {
        *altuena = p3;
        if (p1->y > p2->y) {
            *erdikoa = p1;
            *txikiena = p2;
        } else {
            *erdikoa = p2;
            *txikiena = p1;
        }
    }


}
   // printf("Low point: (%lf, %lf)\n", (*txikiena)->x, (*txikiena)->y);
    //printf("Middle point: (%lf, %lf)\n", (*erdikoa)->x, (*erdikoa)->y);
    //printf("High point: (%lf, %lf)\n", (*altuena)->x, (*altuena)->y);


}


// TODO (transform)
// objektua munduan kokatzen duen matrizetik abiatuta objktuaren erreferentzi sistemara pasatzen duen matrizea lortu
void obtener_CSR_partiendo_de_M(GLdouble* M, GLdouble* MCSR) {
    MCSR[0] = M[0]; MCSR[4] = M[1]; MCSR[8] = M[2]; MCSR [12] = 0;
    MCSR[1] = M[4]; MCSR[5] = M[5]; MCSR[9] = M[6]; MCSR [13] = 0;
    MCSR[2] = M[8]; MCSR[6] = M[9]; MCSR[10] = M[10]; MCSR [14] = 0;
    MCSR[3] = 0;    MCSR[7] = 0;    MCSR[11] = 0;     MCSR [15] = 1;
}

 void kalkulatuH(double H[3],double V[3] , double L[3]){

    H[0] = V[0] + L[0];
    H[1] = V[1] + L[1];
    H[2] = V[2] + L[2];
    normalize(H);
 }
int ikusgaif(double L_[3], double F[3], double angelua){

//printf("cos -> %f\n",dot_product(F, L_));
//printf("cos(2) -> %f\n",(angelua));

if(dot_product(F, L_)>(cos(angelua))){
    return 1;
} else{
return 0;
}
}


void dibujar_triangulo(triobj *optr, int ti, double modelview[16])
{
hiruki *tptr;

punto *goiko, *beheko, *erdiko;
float x1,y1,z1,u1,v1,x2,y2,z2,u2,v2,x3,y3,z3,u3,v3;
float c1x,c1z,c1u,c1v,c2x,c2z,c2u,c2v;
float y; // \in  (-1,1)
float lerrotartea,cambio1,cambio1z,cambio1u,cambio1v,cambio2,cambio2z,cambio2u,cambio2v;
int lerrokop, render;
punto p1,p2,p3,p4,p5,p6,p7;

//p1,p2,p3 munduko erref, p4 normala, p5p6p7 proiekzioa

//printf("hirukian\n");

if (ti >= optr->num_triangles) return;
tptr = optr->triptr+ti;

// printf("p1 puntuak errferentzian: (%.4f, %.4f, %.4f)\n", tptr->p1.x, tptr->p1.y, tptr->p1.z);
// printf("p2 puntuak errferentzian: (%.4f, %.4f, %.4f)\n", tptr->p2.x, tptr->p2.y, tptr->p2.z);
// printf("p3 puntuak errferentzian: (%.4f, %.4f, %.4f)\n", tptr->p3.x, tptr->p3.y, tptr->p3.z);
  //  printMatrix("Mesa", Mesa);
  //  printMatrix("biderkatzea->modelview", modelview);

mxp(&p1,modelview,tptr->p1);
mxp(&p2,modelview,tptr->p2);
mxp(&p3,modelview,tptr->p3);

// printf("p1 puntuak errferentzian: (%.4f, %.4f, %.4f)\n", p1.x, p1.y, p1.z);
// printf("p2 puntuak errferentzian: (%.4f, %.4f, %.4f)\n", p2.x, p2.y, p2.z);
// printf("p3 puntuak errferentzian: (%.4f, %.4f, %.4f)\n", p3.x, p3.y, p3.z);

 if ((p1.z > -0.2|| p2.z > -0.2 || p3.z > -0.2) && perspektiba==1) {//atzetik
 printf("hurbilegi/atzetik\n");
 hurbilegi=1;
     return;
 }
hurbilegi=0;
//printMatrix("biderkatzea Mp", proiekMatrizea);

p4.x=tptr->p1.x + tptr->N[0];
p4.y=tptr->p1.y + tptr->N[1];
p4.z=tptr->p1.z + tptr->N[2];

mxp(&p5, proiekMatrizea, p1); 
mxp(&p6, proiekMatrizea, p2);
mxp(&p7, proiekMatrizea, p3);

// printf("p1 puntuak errferentzian: (%.4f, %.4f, %.4f)\n", p5.x, p5.y, p5.z);
// printf("p2 puntuak errferentzian: (%.4f, %.4f, %.4f)\n", p6.x, p6.y, p6.z);
// printf("p3 puntuak errferentzian: (%.4f, %.4f, %.4f)\n", p7.x, p7.y, p7.z);

double normalVector[3];
double normalaKameraErref[3];
double normalaKameraErrefMp[3];

normalVector[0]=tptr->N[0];
normalVector[1]=tptr->N[1];
normalVector[2]=tptr->N[2];

multMatrixVector3(normalaKameraErref, modelview, normalVector);
multMatrixVector3(normalaKameraErrefMp, proiekMatrizea, normalaKameraErref);



        //     printf("Hau triangeluan->%u,%u,%u \n",optr->color[0],optr->color[1], optr->color[2]);



        if(perspektiba == 0){         
//normala x modelview kamera erreferentzian izateko-----> Datuak beti erreferentzi sistema berdina!!
                                    if(normalaKameraErrefMp[2]<0){//triangulo rojo
                                    glColor3f(1.0, 0.0, 0.0); // Rojo
                                    render=0;
                                    }else{
                                    glColor3f(0.0, 1.0, 0.0); // Verde
                                    render=1;                                              }
            }else{ // z aldrebez?


                                     double cos_theta = cos_angle_camera_normal(modelview, tptr); // Hemen barruan normala erreferentzia sistemaz aldatzeu. Ezta erabiltzen lehenagokoa.
                                        if(cos_theta < 0) {                                   
                                        glColor3f(1.0, 0.0, 0.0); // Rojo
                                        render=0;
                                        }else{
                                        glColor3f(0.0, 1.0, 0.0); // Verde
                                        render=1;
                                        }

        }

        if (backculling == 1 && render==0) {
         // printf("Konkretuki, hau ez da marraztuko\n");
        return;
        }
        else if(backculling == 0 && render==0){


             glBegin(GL_POLYGON);
            glVertex3d(p5.x, p5.y, p5.z);
            glVertex3d(p6.x, p6.y, p6.z);
            glVertex3d(p7.x, p7.y, p7.z);
            glEnd();
            
                if(normala == 1){
                glBegin(GL_LINES);
                glVertex3d(p5.x, p5.y, p5.z);
                glVertex3d(p4.x, p4.y, p4.z);
                glEnd();
                 } 
        return;
         } 

         if(normala == 1){
            
            glBegin(GL_LINES);
            glVertex3d(p5.x, p5.y, p5.z);
            glVertex3d(p4.x, p4.y, p4.z);
            glEnd();
         }  
        if (lineak == 1)
        {
            glBegin(GL_POLYGON);
            glVertex3d(p5.x, p5.y, p5.z);
            glVertex3d(p6.x, p6.y, p6.z);
            glVertex3d(p7.x, p7.y, p7.z);
            glEnd();
             return;
        }
        double Iar=argiGlobala*Ka[0];
        double Iag=argiGlobala*Ka[1];
        double Iab=argiGlobala*Ka[2];

        double Ier;
        double Ieg;
        double Ieb;

        double Ibr;
        double Ibg;
        double Ibb; 

        double Ior;
        double Iog;
        double Iob; 

        double Ikr;
        double Ikg;
        double Ikb;       
     

    if(EguzkiaPiztuta){
       double L[3];
       double angelu;
        multMatrixVector3(L,Mesa, eguzkia->eguzkibektore);
        //printf("%f,%f,%f",L[0], L[1], L[2] );
        normalize(L);
        if (angelu < 0) {
            angelu = 0;
        }
        angelu=dot_product(L,normalaKameraErrefMp);
        if(angelu<0) angelu=0;
        //p5->Mp*Mv*p1
        double V[3]={
            -p5.x,
            -p5.y,
            -p5.z
        };
        double H[3];
         kalkulatuH(H,V,L);
        double angelu2=dot_product(H,normalaKameraErrefMp);
        double lag=pow((angelu2 ), distira);// ondo...?

         Ier=((angelu)*Kd[0] + lag*Ks[0])*eguzkia->intentsitate[0];
         Ieg=((angelu)*Kd[1] + lag*Ks[1])*eguzkia->intentsitate[1];
         Ieb=((angelu)*Kd[2] + lag*Ks[2])*eguzkia->intentsitate[2];

    }
    if(BonbilaPiztuta){        //p5->Mp*Mv*p1

            double L[3];
            double bonbilaKam[3];
            double angelu;
        multMatrixVector3(bonbilaKam,Mesa,bonbila->puntua);
        L[0]= bonbilaKam[0]-p5.x;
        L[1]=bonbilaKam[1]-p5.y;
        L[2]=bonbilaKam[2]-p5.z;
        angelu=dot_product(L,normalaKameraErrefMp);
        if(angelu<0) angelu=0;
        double V[3]={
            -p5.x,
            -p5.y,
            -p5.z
        };
        double H[3];
         kalkulatuH(H,V,L);
        double angelu2=dot_product(H,normalaKameraErrefMp);
        double lag=pow((angelu2), distira);
         Ibr=((angelu)*Kd[0] + lag*Ks[0])*bonbila->intentsitate[0];
         Ibg=((angelu)*Kd[1] + lag*Ks[1])*bonbila->intentsitate[1];
         Ibb=((angelu)*Kd[2] + lag*Ks[2])*bonbila->intentsitate[2];

        // printf("Intentsitateak Ibonbila totala(R,G,B)->%f,%f,%f \n", Ibr, Ibg,Ibb);

    }
    if(ObjektufokuaPiztuta){ //p5->Mp*Mv*p1
        double L[3];
        double L_[3];
    
        double modelviewLFokua[16];

        double fokua[3];
        double fokuaMv[3];
        double fokuaMvMp[3];

        double fokuaDir[3];
        double fokuaDirMv[3];
        double F[3];

        fokua[0]= sel_ptr->triptr->p1.x;
        fokua[1]= sel_ptr->triptr->p1.y;
        fokua[2]= sel_ptr->triptr->p1.z;

        fokuaDir[0]=0;
        fokuaDir[1]=0;
        fokuaDir[2]=1;

        matrizeaKalk(modelviewLFokua, Mesa, sel_ptr->mptr->m); // FOKUAren  kamararen erreferentzia sistema pasatzeko matrizea

        multMatrixVector3(fokuaMv, modelviewLFokua, fokua);
        multMatrixVector3(fokuaMvMp,proiekMatrizea, fokuaMv);

        multMatrixVector3(fokuaDirMv, modelviewLFokua, fokuaDir);
        multMatrixVector3(F,proiekMatrizea, fokuaDirMv);

        normalize(F);        


        L[0]=fokuaMvMp[0]-p5.x;
        L[1]=fokuaMvMp[1]-p5.y;
        L[2]=fokuaMvMp[2]-p5.z;

        normalize(L);        

        
        L_[0]=-L[0];
        L_[1]=-L[1];
        L_[2]=-L[2];

        int ikusgai=ikusgaif(L_,F,angeluaFokuaObjektua);
        //ikusgai=1;
       // printf("Ikusgai->%d",ikusgai);
        double angelu=dot_product(L,normalaKameraErrefMp);
        if(angelu<0) angelu=0;
        double V[3]={
            -p5.x,
            -p5.y,
            -p5.z
        };
        double H[3];
         kalkulatuH(H,V,L);
        double angelu2=dot_product(H,normalaKameraErrefMp);
        double lag=pow((angelu2), distira);
         Ior=((angelu)*Kd[0]+ lag*Ks[0])*IObjFokua*ikusgai;
         Iog=((angelu)*Kd[1] + lag*Ks[1])*IObjFokua*ikusgai;
         Iob=((angelu)*Kd[2] + lag*Ks[2])*IObjFokua*ikusgai;
        //printf("Intentsitateak Iobj totala(R,G,B)->%f,%f,%f \n", Ior, Iog,Iob);

    }if(KamerafokuaPiztuta){
        double F[3];
        double L[3];
        double L_[3];
        double angelu;

        F[0]=0;
        F[1]=0;
        F[2]=-1;

        L[0]=-p5.x;
        L[1]=-p5.y;
        L[2]=-p5.z;

        L_[0]=-L[0];
        L_[1]=-L[1];
        L_[2]=-L[2];

        int ikusgai=ikusgaif(L_,F,angeluaFokuaKamera);

        angelu=dot_product(L,normalaKameraErrefMp);
        if(angelu<0) angelu=0;
         double V[3]={
            -p5.x,
            -p5.y,
            -p5.z
        };
        double H[3];
        kalkulatuH(H,V,L);
        double angelu2=dot_product(H,normalaKameraErrefMp);
        double lag=pow((angelu2), distira);
         Ikr=((angelu)*Kd[0] + lag*Ks[0])*IKamFokua*ikusgai;
         Ikg=((angelu)*Kd[1] + lag*Ks[1])*IKamFokua*ikusgai;
         Ikb=((angelu)*Kd[2] + lag*Ks[2])*IKamFokua*ikusgai;

    }
    //printf("Intentsitateak Ieguzkia totala(R,G,B)->%f,%f,%f \n", Ibr, Ibg,Ibb);

    double Ir=Iar + Ier + Ibr + Ior + Ikr;
    double Ig=Iag + Ieg + Ibg + Iog + Ikg;
    double Ib=Iab + Ieb + Ibb + Iob + Ikb;

    if(Ir>1.0) Ir=1.0;
    if(Ig>1.0) Ig=1.0; 
    if(Ib>1.0) Ib=1.0;
     //printf("Intentsitateak I totala(R,G,B)->%f,%f,%f \n", Ir, Ig,Ib);

    //printf("Intentsitateak tptr(R,G,B)->%d,%d,%d \n",  tptr->color[0],  tptr->color[1], tptr->color[2]);

                 tptr->color[0]= Ir*optr->color[0];
                 tptr->color[1]= Ig*optr->color[1];
                 tptr->color[2]= Ib*optr->color[2];

  // printf("Intentsitateak tptr gero(R,G,B)->%d,%d,%d \n",  tptr->color[0],  tptr->color[1], tptr->color[2]);

lerrotartea = 2.0/(float)dimentsioa;
ordenatu(&p5,&p6,&p7, &beheko, &erdiko, &goiko);
//printf("goiko %.2f\n", goiko->y);
//printf("erdiko %.2f\n", erdiko->y);
//printf("beheko %.2f\n", beheko->y);
     if (goiko->y == beheko->y) {
        // Altura berean daude
    dibujar_linea_z(erdiko->y, beheko->x, beheko->z, beheko->u, beheko->v, goiko->x, goiko->z, goiko->u, goiko->v,tptr->color);
         return;

    }
     for (y = beheko->y; y <= goiko->y; y += lerrotartea) {
        if (y <= erdiko->y) {
            //behetik erdira
            c1x = beheko->x + (y - beheko->y) / (erdiko->y - beheko->y) * (erdiko->x - beheko->x);
            c1u = beheko->u + (y - beheko->y) / (erdiko->y - beheko->y) * (erdiko->u - beheko->u);
            c1v = beheko->v + (y - beheko->y) / (erdiko->y - beheko->y) * (erdiko->v - beheko->v);
            c1z = beheko->z + (y - beheko->y) / (erdiko->y - beheko->y) * (erdiko->z- beheko->z);


            c2x = beheko->x + (y - beheko->y) / (goiko->y - beheko->y) * (goiko->x - beheko->x);
            c2u = beheko->u + (y - beheko->y) / (goiko->y - beheko->y) * (goiko->u - beheko->u);
            c2v = beheko->v + (y - beheko->y) / (goiko->y - beheko->y) * (goiko->v - beheko->v);
            c2z = beheko->z + (y - beheko->y) / (goiko->y - beheko->y) * (goiko->z- beheko->z);

        } else {
            // erditik gora
            c1x = erdiko->x + (y - erdiko->y) / (goiko->y - erdiko->y) * (goiko->x - erdiko->x);
            c1u = erdiko->u + (y - erdiko->y) / (goiko->y - erdiko->y) * (goiko->u - erdiko->u);
            c1v = erdiko->v + (y - erdiko->y) / (goiko->y - erdiko->y) * (goiko->v - erdiko->v);
            c1z = erdiko->z + (y - erdiko->y) / (goiko->y - erdiko->y) * (goiko->z- erdiko->z);


            c2x = beheko->x + (y - beheko->y) / (goiko->y - beheko->y) * (goiko->x - beheko->x);
            c2u = beheko->u + (y - beheko->y) / (goiko->y - beheko->y) * (goiko->u - beheko->u);
            c2v = beheko->v + (y - beheko->y) / (goiko->y - beheko->y) * (goiko->v - beheko->v);
            c2z = beheko->z + (y - beheko->y) / (goiko->y - beheko->y) * (goiko->z- beheko->z);

        }
        //ondo marrazteko beharrezkoa, c1x puntua ezkerrean egon behar da.
 	 if (c1x > c2x) {
            float temp_x = c1x, temp_u = c1u, temp_v = c1v;
            c1x = c2x; c1u = c2u; c1v = c2v;
            c2x = temp_x; c2u = temp_u; c2v = temp_v;
        }
        // y balioan marraztu
        dibujar_linea_z(y, c1x, c1z, c1u, c1v, c2x, c2z, c2u, c2v,tptr->color);
    }
 }


static void marraztu(void) {
    int i;
    triobj *auxptr;


    // Verificar que hay objetos para dibujar
    if (foptr == 0) return;

    // Limpiar el buffer de color y profundidad según las condiciones
    if (objektuak == 1) {
        glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
    } else {
        if (denak == 0) {
        glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
        }
    }
   // printMatrix("Kamararen matrizea", kamptr->mptr->m);

    MesaLortu();  // mesa aldagai globala da


    mpkalkulatu(); //proiekMatrizea aldagai globala izando da


    printf("\n");
        printf("\n");



 printf("Perspektiba:");
printf("%d", perspektiba);
        printf("\n");
if(argia==0){
if(kamera == 1){
    printf("Kamera aldatzen...:\n");

}
else{
    printf("Objektua aldatzen...\n");

}}
else{
    printf("Argia aldatzen:\n");
    
        switch (zeinArgi) {
        case 0:
            printf("Eguzkia...\n");
            break;
        case 1:
            printf("Bonbila...\n");
            break;
        case 2:
            printf("Objektuan fokua...\n");
            break;
        case 3:
            printf("Kameran fokua...\n");
            break;

               }

}

 printf("AldaketaLokala:");
printf("%d", ald_lokala);
printf("\n");

 printf("Hegaldi modua:");
printf("%d", hegaldi);
printf("\n");

printf("Backculling:");
printf("%d\n", backculling);

printf("Objektu ikuspegi:");
printf("%d\n", objektuIkuspegi);

printf("Normala:");
printf("%d\n", normala);

printf("Eguzkia:");
printf("%d\n", EguzkiaPiztuta);

printf("Bonbila:");
printf("%d\n", BonbilaPiztuta);

printf("Kamararen argia:");
printf("%d\n", KamerafokuaPiztuta);

printf("Objektuaren argia:");
printf("%d\n", ObjektufokuaPiztuta);

    if (objektuak == 1) {   
        if (denak == 1) {
            for (auxptr = foptr; auxptr != 0; auxptr = auxptr->hptr) {

                double modelview[16];
               matrizeaKalk(modelview, Mesa, auxptr->mptr->m); 

                for (i = 0; i < auxptr->num_triangles; i++) {
                    dibujar_triangulo(auxptr, i,modelview);
                }
            }
        } else {
            auxptr = sel_ptr;
            if (auxptr != 0) {
                double modelview[16];
                matrizeaKalk(modelview, Mesa, auxptr->mptr->m);

                for (i = 0; i < auxptr->num_triangles; i++) {
                    dibujar_triangulo(auxptr, i,modelview);
                }
            }
        }
    } else {

        auxptr = sel_ptr;
        if (auxptr != 0) {
            double modelview[16];
            matrizeaKalk(modelview, Mesa, auxptr->mptr->m);

                    dibujar_triangulo(auxptr, i,modelview);
        }
    }
    
    glFlush(); 
}



void x_aldaketa(int dir) {
    double angulo = 3.0 * dir;
    double rad = angulo * 3.1415 / 180.0;

    double rotazio_x[16] = {
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1
    };

if(argia==1 && zeinArgi==0){

    rotazio_x[5] = cos(rad);
        rotazio_x[6] = sin(rad);
        rotazio_x[9] = -sin(rad);
        rotazio_x[10] = cos(rad);

        objektuari_aldaketa_sartu_esk(rotazio_x); // Transformazio lokala

}


   if(kamera == 1 && hegaldi ==1){
    
 // rotazioa
        rotazio_x[5] = cos(rad);
        rotazio_x[6] = sin(rad);
        rotazio_x[9] = -sin(rad);
        rotazio_x[10] = cos(rad);

        objektuari_aldaketa_sartu_esk(rotazio_x); // Transformazio lokala

    }else{ 

    if (aldaketa == 'r') {
        // rotazioa
        rotazio_x[5] = cos(rad);
        rotazio_x[6] = sin(rad);
        rotazio_x[9] = -sin(rad);
        rotazio_x[10] = cos(rad);
    } else {
        //Translazioa 
        rotazio_x[12] += dir * 0.1;  
    }

    if (ald_lokala) {
        objektuari_aldaketa_sartu_esk(rotazio_x); // lokala

    } else {
        objektuari_aldaketa_sartu_ezk(rotazio_x); //  globala

    }
}
}

void y_aldaketa(int dir) {
    double angulo = 3.0 * dir;
    double rad = angulo * 3.1415 / 180.0;

    double rotazio_y[16] = {
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1
    };
if(argia==1 && zeinArgi==0){

        rotazio_y[0] = cos(rad);
        rotazio_y[2] = -sin(rad);
        rotazio_y[8] = sin(rad);
        rotazio_y[10] = cos(rad);
   objektuari_aldaketa_sartu_esk(rotazio_y);
}else{


   if(kamera == 1 && hegaldi == 1){


        rotazio_y[0] = cos(rad);
        rotazio_y[2] = -sin(rad);
        rotazio_y[8] = sin(rad);
        rotazio_y[10] = cos(rad);
   objektuari_aldaketa_sartu_esk(rotazio_y);

    }else{

    if (aldaketa == 'r') {
        rotazio_y[0] = cos(rad);
        rotazio_y[2] = -sin(rad);
        rotazio_y[8] = sin(rad);
        rotazio_y[10] = cos(rad);
    } else {
        rotazio_y[13] += dir * 0.1;  
    }

    if (ald_lokala) {
        objektuari_aldaketa_sartu_esk(rotazio_y);

    } else {
        objektuari_aldaketa_sartu_ezk(rotazio_y);

    }
    }
}
}   

void z_aldaketa(int dir) {
    double angulo = 3.0 * dir;
    double rad = angulo * 3.1415 / 180.0;

    double rotazio_z[16] = {
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1
    };
    if(argia==1){
                rotazio_z[0] = cos(rad);
        rotazio_z[1] = sin(rad);
        rotazio_z[4] = -sin(rad);
        rotazio_z[5] = cos(rad);
        objektuari_aldaketa_sartu_esk(rotazio_z); 


    }else{
            if (aldaketa == 'r') {
        rotazio_z[0] = cos(rad);
        rotazio_z[1] = sin(rad);
        rotazio_z[4] = -sin(rad);
        rotazio_z[5] = cos(rad);
    } else {
        rotazio_z[14] = dir * 0.1; 
    }

    if (ald_lokala) {
        objektuari_aldaketa_sartu_esk(rotazio_z); 

    } else {
        objektuari_aldaketa_sartu_ezk(rotazio_z); 

    }

    }


}


void kameraMugituZ(int dir) { // 
    // Calcular el cambio en el eje Z
    double angulo = 3.0 * dir;
    double rad = angulo * 3.1415 / 180.0;
    
   
    double rotazio_z[16] = {
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1
    };
    rotazio_z[14] = dir * 0.1;

    objektuari_aldaketa_sartu_esk(rotazio_z);

    if(hegaldi == 0){
    lookAt(sel_ptr->mptr->m[12], sel_ptr->mptr->m[13], sel_ptr->mptr->m[14]);
    }

}




void undo() { // HPTR=0 

    if(kamera == 0){
        if(sel_ptr->mptr->hptr==0){
        printf("Ezinezkoa(objektua)!!");
        fflush(stdout); //bufferra hustu
        return;
    }
        if (sel_ptr && sel_ptr->mptr && sel_ptr->mptr->hptr) {
            // aurreko matrizeari pointerra jarri
            mlist *aurreko = sel_ptr->mptr->hptr;

            //oraingo matrizea kendu
            free(sel_ptr->mptr);

            // aurreko matrizea berreskuratu
            sel_ptr->mptr = aurreko;
        }

    }else{
         if(kamptr->mptr->hptr==0){
        printf("Ezinezkoa(kamera)!!");
        fflush(stdout); //bufferra hustu
        return;
    }
        if (kamptr && kamptr->mptr && kamptr->mptr->hptr) {
            // aurreko matrizeari pointerra jarri
            mlist *aurreko = kamptr->mptr->hptr;

            //oraingo matrizea kendu
            free(kamptr->mptr);

            // aurreko matrizea berreskuratu
            kamptr->mptr = aurreko;
        }

        
    }

}


void read_from_file(char *fitx)
{
    int i,retval;
    triobj *optr;
    printf("%s fitxategitik datuak hartzera\n",fitx);
    optr = (triobj *)malloc(sizeof(triobj));
    unsigned char *rgbptr = NULL;
    //TODO (transform...)
    //retval = cargar_triangulos(fitx, &(optr->num_triangles), &(optr->triptr));
    retval = cargar_triangulos_color(fitx, &(optr->num_triangles), &(optr->triptr), &rgbptr);
    printf("Hemen emaitza retval:%d\n",retval );

    if (retval ==-1)
         {
         printf("%s fitxategitik datuak hartzerakoan arazoak izan ditut\n",fitxiz);
         free(optr);
         }
       else if(retval==15 || retval==9)
         {
         triangulosptr = optr->triptr;
         printf("objektuaren matrizea...\n");
         optr->mptr = (mlist *)malloc(sizeof(mlist));
         for (i=0; i<16; i++) optr->mptr->m[i] =0;
         optr->mptr->m[0] = 1.0;
         optr->mptr->m[5] = 1.0;
         optr->mptr->m[10] = 1.0;
         optr->mptr->m[15] = 1.0;
         optr->mptr->hptr = 0; //hurrengoa da 0
         printf("objektu zerrendara doa informazioa...\n");
         optr->hptr = foptr;
         foptr = optr;
         sel_ptr = optr;

        printf("Kamera hasieratzen...\n");
        kamptr = (triobj *)malloc(sizeof(triobj));  
        kamptr->mptr = (mlist *)malloc(sizeof(mlist));
         for (i=0; i<16; i++) kamptr->mptr->m[i] =0;
         kamptr->mptr->m[0] = 1.0;
         kamptr->mptr->m[5] = 1.0;
         kamptr->mptr->m[10] = 1.0;
         kamptr->mptr->m[15] = 1.0;
        kamptr->mptr->m[14] = 2.5;

         kamptr->mptr->hptr = 0; //hurrengoa da 0triobj *kamptr;

        for (i=0; i<16; i++) proiekMatrizea[i] =0;
        proiekMatrizea[0] = 1.0;
        proiekMatrizea[5] = 1.0;
        proiekMatrizea[10] = 1.0;
        proiekMatrizea[15] = 1.0;

         for (i=0; i<16; i++) Mesa[i] =0;
        Mesa[0] = 1.0;
         Mesa[5] = 1.0;
         Mesa[10] = 1.0;
         Mesa[15] = 1.0;

         
        for (i = 0; i < optr->num_triangles; i++) {
            bektoreNormalaKalk(&(optr->triptr[i]));
        }
         printf("Argiak hasieratzen...\n");
        eguzkia = (eguzkiobj *)malloc(sizeof(eguzkiobj));  

        for (i=0; i<3; i++){
             eguzkia->intentsitate[i]=0;
             eguzkia->eguzkibektore[i]=0;
        } 
        eguzkia->eguzkibektore[1]=1;
        eguzkia->intentsitate[0]=1;

        bonbila = (bonbilaobj *)malloc(sizeof(bonbilaobj));  

        for (i=0; i<3; i++) {
            bonbila->intentsitate[i]=0.0;
            bonbila->puntua[i]=0;

        }
        bonbila->puntua[0]=-0.5;
        bonbila->puntua[1]=0.3;
        bonbila->puntua[2]=0.0;
        bonbila->intentsitate[1]=0.8;
    

        printf("Plastiko distiratsuan konstanteekin egingo dugu lan\n");
        for(int i=0;i<3;i++){
            Ka[i]=0.2;
            Kd[i]=0.8;
            Ks[i]=0.5;
        }   

         }
        if (retval == 9 && rgbptr != NULL) {
           // printf("Koloriak: R=%d, G=%d, B=%d\n", rgbptr[0], rgbptr[1], rgbptr[2]);
                optr->color[0] = rgbptr[0] ; 
                optr->color[1] = rgbptr[1] ;
                optr->color[2] = rgbptr[2] ;
            for (i = 0; i < optr->num_triangles; i++) {
                optr->triptr[i].color[0] = rgbptr[0] ; 
                optr->triptr[i].color[1] = rgbptr[1] ;
                optr->triptr[i].color[2] = rgbptr[2] ;
            }
        }
     printf("datuak irakurrita\n");
}

// This function will be called whenever the user pushes one key
static void teklatua (unsigned char key, int x, int y)
{
int retval;
int i;
FILE *obj_file;

switch(key)
	{
	case 13:
	        if (foptr != 0)  // objekturik ez badago ezer ez du egin behar
	                         // si no hay objeto que no haga nada
	            {
	            indexx ++;  // azkena bada lehenengoa bihurtu
		                // pero si es el último? hay que controlarlo!
		    if (indexx == sel_ptr->num_triangles)
		        {
		        indexx = 0;
		        if ((denak == 1) && (objektuak == 0))
		            {
		            glClear( GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT );
		            glFlush();
		            }
		        }
		    }
		break;
	case 'd':
		if (denak == 1) denak = 0;
		    else denak = 1;
		break;
    case 'C':
        if(objektuIkuspegi==0){
    kamera=0;
    objektuIkuspegi=1; 
        }else{
           objektuIkuspegi = 0;
        } 
        break;
	case 'o':
		if (objektuak == 1) objektuak = 0;
		    else objektuak = 1;
		break;
	case 'c':
		if (argia == 1  && objektuIkuspegi == 0){
            kamera = 0;
            hegaldi = 0;
            objektu=1;
            argia=0;
        } 
		else if(objektuIkuspegi == 0 && kamera == 0 && objektu==1 ){
        objektuIkuspegi=0;
        kamera = 1; 
        objektu=0;
        hegaldi= 1;
        argia=0;
        }else if(kamera == 1){
         kamera=0;
         argia=1;
        }
		break;
	case 'l':
		if (lineak == 1) lineak = 0;
		    else lineak = 1;
		break;
	case 't':
	        aldaketa = 't';
        if(kamera == 1){
            printf("Kameran beti dira mugimendu berdinak!");
        }
		break;
	case 'r':
		aldaketa = 'r';
        if(kamera == 1){
            printf("Kameran beti dira mugimendu berdinak!");
        }
		break;
	case 'g':
    if (ald_lokala == 1 && kamera == 0) {
        ald_lokala = 0; 
    } else if (ald_lokala == 0 && kamera == 0) {
        ald_lokala = 1;  
    } else if(kamera ==1 && hegaldi == 0){
        hegaldi =1;
        kamera=1;
         ald_lokala = 1;
    }else{
      hegaldi = 0;
                   lookAt(sel_ptr->mptr->m[12],sel_ptr->mptr->m[13],sel_ptr->mptr->m[14]);

        kamera=1;
         ald_lokala = 1;

    }
		break;
    case 'G':
		if (kamera==1) hegaldi = 1;
		break;
               //p eta b 

case 'x':
if(argia==0){
       if (kamera == 1 && hegaldi == 0) {
        kameraAnalisi(1, 'x');  

    } else {
        x_aldaketa(1); 
    }

}else{
    if(zeinArgi==0){
    x_aldaketa(1);
    }else if(zeinArgi==1){
    bonbila->puntua[0]+=0.1;
    }
    
}
 
    break;

case 'y':
if(argia==0){
        if (kamera == 1 && hegaldi == 0) {
        kameraAnalisi(1, 'y');
    } else {
        y_aldaketa(1); 
    }

}else{

  if(zeinArgi==0){
    y_aldaketa(1);
    }else if(zeinArgi==1){
    bonbila->puntua[1]+=0.1;
    }

}

    break;

case 'z':
if(argia==0){
       if(kamera==0){
        z_aldaketa(1);
        }else if(kamera == 1){
            kameraMugituZ(1);

            }

}else{

      if(zeinArgi==0){
        printf("Ezinezkoa!\n");
    }else if(zeinArgi==1){
    bonbila->puntua[2]+=0.1;
    }

}
    break;


case 'X':

if(argia==0){
           if (kamera == 1 && hegaldi == 0) {
        kameraAnalisi(-1, 'x'); 
    } else {
        x_aldaketa(-1); 
    }
}else{

        if(zeinArgi==0){
    x_aldaketa(1);
    }else if(zeinArgi==1){
      bonbila->puntua[0]-=0.1;
    }

}
    break;
        case 'Y':
        if(argia==0){
              if (kamera == 1 && hegaldi == 0) {
        kameraAnalisi(-1, 'y');  
    } else {
        y_aldaketa(-1); 
    }

        }else{

       if(zeinArgi==0){
    y_aldaketa(1);
    }else if(zeinArgi==1){
    bonbila->puntua[1]-=0.1;
    }


}

    break;        
        case 'Z':
        if(argia==0){
                    if(kamera==0){
        z_aldaketa(-1);
            }else if(kamera == 1){
                kameraMugituZ(-1);

            }

        }else{
            if(zeinArgi==0){
        printf("Ezinezkoa!\n");
          }else if(zeinArgi==1){
            bonbila->puntua[2]-=0.1;
         }

        }
   
               
        break;
        case 'p':
        if(perspektiba == 0){
        perspektiba=1;
        }
            
        else{
    perspektiba = 0;
        }
            break;
        case 'u':
                undo();
                break;

        case '1':
        if(EguzkiaPiztuta==0){
            EguzkiaPiztuta=1;
        }else{
            EguzkiaPiztuta=0;
        }
         break;
         
        case '2':
        if(BonbilaPiztuta==0){
            BonbilaPiztuta=1;
        }else{
            BonbilaPiztuta=0;
        }
         break;

          case '4':
        if(KamerafokuaPiztuta==0){
            KamerafokuaPiztuta=1;
        }else{
            KamerafokuaPiztuta=0;
        }
         break;

          case '3':
        if(ObjektufokuaPiztuta==0){
            ObjektufokuaPiztuta=1;
        }else{
            ObjektufokuaPiztuta=0;
        }
         break;
         
	case 'f':
	        /*Ask for file*/
	        printf("idatzi fitxategi izena\n");
	        scanf("%s", &(fitxiz[0]));
	        read_from_file(fitxiz);
	        indexx = 0;
                break;

    case 'b':
     if(backculling==0){
        backculling=1;
            }else {
        backculling=0;
            }
         break;      
    
       /* case 'S':  // save to file
	        printf("idatzi fitxategi izena\n");
	        scanf("%s", &(fitxiz[0]));
                if ((obj_file = fopen(fitxiz, "w")) == NULL)
                         {
                         printf("ezin fitxategia ireki\n");
                         }
                     else
                         {
                         for (i =0; i < sel_ptr->num_triangles; i++)
                            {
                            fprintf(obj_file,"t %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
                                 sel_ptr->triptr[i].p1.x-250, sel_ptr->triptr[i].p1.y-250, sel_ptr->triptr[i].p1.z,
                                 sel_ptr->triptr[i].p1.u, sel_ptr->triptr[i].p1.v,
                                 sel_ptr->triptr[i].p2.x-250, sel_ptr->triptr[i].p2.y-250, sel_ptr->triptr[i].p2.z,
                                 sel_ptr->triptr[i].p2.u, sel_ptr->triptr[i].p2.v,
                                 sel_ptr->triptr[i].p3.x-250, sel_ptr->triptr[i].p3.y-250, sel_ptr->triptr[i].p3.z,
                                 sel_ptr->triptr[i].p3.u, sel_ptr->triptr[i].p3.v );
                            }
                         fclose(obj_file);
                         }
                break; */
       case 9: /* <TAB> para cambiar el objeto seleccionado */
       if(argia==0){
            if (foptr != 0) { // Si no hay objetos, no hace nada
        sel_ptr = sel_ptr->hptr; // Avanza al siguiente objeto en la lista
        /* Selección circular: si llegamos al final de la lista, volvemos al primer elemento */
        if (sel_ptr == 0) sel_ptr = foptr;
        indexx = 0; // Reiniciamos el índice del polígono seleccionado

   
    }
     if(hegaldi == 0 && kamera == 1){
        lookAt(sel_ptr->mptr->m[12],sel_ptr->mptr->m[13],sel_ptr->mptr->m[14]);
    }
       }else{
        zeinArgi++; 
        if (zeinArgi == 4) zeinArgi = 0; 
           
       }
    
    break;
    case 127: /* <Supr> para eliminar el objeto seleccionado */
    if (sel_ptr != 0) { // Si hay un objeto seleccionado
        triobj *to_delete = sel_ptr;

        // Si el objeto a eliminar es el primero en la lista
        if (to_delete == foptr) {
            foptr = foptr->hptr; // Actualiza el primer objeto de la lista
            sel_ptr = foptr;      // Nuevo objeto seleccionado
        } else {
            // Recorrer la lista para encontrar el objeto anterior al que queremos eliminar
            triobj *prev = foptr;
            while (prev->hptr != to_delete) {
                prev = prev->hptr;
            }
            prev->hptr = to_delete->hptr; // Enlaza el anterior con el siguiente
            sel_ptr = prev->hptr ? prev->hptr : prev; // Actualizar selección
        }

        // Liberar memoria del objeto eliminado
        free(to_delete->triptr);
        free(to_delete->mptr);
        free(to_delete);

        // Si la lista queda vacía, restablece `sel_ptr` a NULL
        if (foptr == 0) {
            sel_ptr = 0;
            printf("Azken objetua borratu da, beraz, programa bukatu egingo da.\n");
            exit(0);
        } else {
            printf("Objektua borratu da eta memoria askatu da.\n");
           }
        // Limpiar el búfer y refrescar la pantalla
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glutPostRedisplay();  // Llamar a redibujar la pantalla
    }
    break;

    case 'n':
    if(normala == 1){
        normala=0;
    }else{
        normala=1;
    }


    break;

	case 27:  // <ESC>
		exit( 0 );
		break;
        case '+':  
if(argia == 1){

    if(zeinArgi==2){
        angeluaFokuaObjektua+=0.05;
        printf("AngeluaObjektuarena:%f\n",angeluaFokuaObjektua);

    }
    if(zeinArgi==3){
        angeluaFokuaKamera+=0.05;
        printf("Angelua Kamera:%f\n",angeluaFokuaKamera);

    }

}else{
            kamera=0;
 
    eskala_aldaketa(1.1); //10%

}

    break;
case '-':  
if(argia == 1){

    if(zeinArgi==2){
        if(angeluaFokuaObjektua>0.1){
        angeluaFokuaObjektua-=0.05;
        }
        printf("AngeluaObjektuarena:%f\n",angeluaFokuaObjektua);
    }
    if(zeinArgi==3){
        if(angeluaFokuaKamera>0.1){
        angeluaFokuaKamera-=0.05;
        }
        printf("Angelua Kamera:%f\n",angeluaFokuaKamera);

    }

}else{
    kamera=0;
    eskala_aldaketa(0.9); // -10%

}

    break;
    
	default:
		printf("%d %c\n", key, key );
	}

// The screen must be drawn to show the new triangle
glutPostRedisplay();
}


void viewportberria (int zabal, int garai)
{
if (zabal < garai)  dimentsioa = zabal;
    else  dimentsioa = garai;
glViewport(0,0,dimentsioa,dimentsioa);
printf("linea kopuru berria = %d\n",dimentsioa);
}

int main(int argc, char** argv)
{
int retval;

	printf(" Triangeluak: barneko puntuak eta testura\n Triángulos con puntos internos y textura \n");
	printf("Press <ESC> to finish\n");
	glutInit(&argc,argv);
	glutInitDisplayMode ( GLUT_RGB|GLUT_DEPTH );
	dimentsioa = 500;
	glutInitWindowSize ( dimentsioa, dimentsioa );
	glutInitWindowPosition ( 800, 100 );
	glutCreateWindow( "KBG/GO praktika" );

	glutDisplayFunc( marraztu );
	glutKeyboardFunc( teklatua );
	glutReshapeFunc( viewportberria);
	/* we put the information of the texture in the buffer pointed by bufferra. The dimensions of the texture are loaded into dimx and dimy */
        retval = load_ppm("proba.ppm", &bufferra, &dimx, &dimy);
        if (retval ==-1)
            {
            printf("Ez dago testuraren fitxategia (testura.ppm)\n");
            exit(-1);
            }

	glClearColor( 0.0f, 0.0f, 0.7f, 1.0f );
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glEnable(GL_DEPTH_TEST); // activar el test de profundidad (Z-buffer)
        glDepthFunc(GL_GREATER); // Handiena marraztu
        glClearDepth(0.0); // Buferra hasieratzeko balioa 0 izan dadila 1 izan beharrean.
        denak = 1;
        lineak =0;
        objektuak = 1;
        kamera = 0;
        hegaldi = 1;
        foptr = 0;
        sel_ptr = 0;
        aldaketa = 't';
        ald_lokala =1;
        perspektiba =1;
        backculling=1;
        normala=0;
        hurbilegi=0;
        argia=0;
        EguzkiaPiztuta=1;
        BonbilaPiztuta=0;
        ObjektufokuaPiztuta=0;
        KamerafokuaPiztuta=0;
        argiGlobala=0.8; 
        objektu=1;     
        distira=32;
        angeluaFokuaKamera=0.45;
        angeluaFokuaObjektua=0.40;
        IKamFokua=0.7;
        IObjFokua=0.7;
        zeinArgi=0  ;

        if (argc>1) read_from_file(argv[1]);
            else{
                read_from_file("abioia-1+1.txt");
            if (sel_ptr != 0) 
                { sel_ptr->mptr->m[12] = -1.0;
                if (sel_ptr->color !=0)sel_ptr->color[0]=0;
                }   
            read_from_file("abioia-1+1.txt");
            if (sel_ptr != 0) 
                { sel_ptr->mptr->m[12] = 1.0;
                if (sel_ptr->color !=0) sel_ptr->color[1]=0;
                }   
            read_from_file("abioia-1+1.txt");
            if (sel_ptr != 0) 
                { 
                sel_ptr->mptr->m[13] = -0.4;
                sel_ptr->mptr->m[14] = 0.5;
                if (sel_ptr->color !=0) sel_ptr->color[2]=0;
                }   
            read_from_file("abioia-1+1.txt");
            if (sel_ptr != 0) 
                { sel_ptr->mptr->m[14] = -0.7;
                if (sel_ptr->color !=0) sel_ptr->color[2]=0;
                }        
            read_from_file("abioia-1+1.txt");
            if (sel_ptr != 0) 
                { sel_ptr->mptr->m[13] = 0.7;
                if (sel_ptr->color !=0) sel_ptr->color[2]=0;
                }                        
            read_from_file("abioia-1+1.txt");
            } 
	glutMainLoop();

	return 0;
}
     