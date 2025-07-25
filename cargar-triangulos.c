#include <stdio.h>
#include <malloc.h>
#define MAXLINE 200

/*
typedef struct punto
{
float x, y, z, u,v;
} punto;

typedef struct hiruki
{
punto p1,p2,p3;
} hiruki;
*/
#include "cargar-triangulo.h"

int cargar_triangulos(char *fitxiz, int *hkopptr, hiruki **hptrptr)
{
 FILE *obj_file;
 char line[MAXLINE];
 int i, num_triangles;

 if ((obj_file = fopen(fitxiz, "r")) == NULL) 
        {
        *hkopptr= 0;
        return(-1);
        }
 num_triangles=0;
 while (fscanf(obj_file, "\n%[^\n]", line) > 0) 
	{
        if (line[0] == 't')  // triangulo!
		{
                 num_triangles++;
		}
	}
 fclose(obj_file);
 *hkopptr= num_triangles;
 *hptrptr = (hiruki *)malloc(num_triangles * sizeof (hiruki));
 
 obj_file = fopen(fitxiz, "r");

 i=0;
 while (fscanf(obj_file, "\n%[^\n]", line) > 0) 
	{
        if (line[0] == 't')  // triangulo!
		{
		sscanf(line + 2, "%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f", &((*hptrptr)[i].p1.x),&((*hptrptr)[i].p1.y),&((*hptrptr)[i].p1.z),
		                                                   &((*hptrptr)[i].p1.u),&((*hptrptr)[i].p1.v),
		                                                   &((*hptrptr)[i].p2.x),&((*hptrptr)[i].p2.y),&((*hptrptr)[i].p2.z),
		                                                   &((*hptrptr)[i].p2.u),&((*hptrptr)[i].p2.v),
		                                                   &((*hptrptr)[i].p3.x),&((*hptrptr)[i].p3.y),&((*hptrptr)[i].p3.z),
		                                                   &((*hptrptr)[i].p3.u),&((*hptrptr)[i].p3.v));
                i++;
		}
	}
 fclose(obj_file);
 return(1);
}


int cargar_triangulos_color(char *fitxiz, int *hkopptr, hiruki **hptrptr, unsigned char **rgbptr) {
    FILE *obj_file;
    char line[MAXLINE];
    int i, num_triangles;
    int zkop;
    int r = 255, g = 255, b = 255, kolorea = 0;
    float v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15;

    if ((obj_file = fopen(fitxiz, "r")) == NULL) {
        *hkopptr = 0;
        printf("Error: no se pudo abrir el archivo %s\n", fitxiz);
        return -1;
    }

    // Contar el número de triángulos
    num_triangles = 0;
    while (fscanf(obj_file, "\n%[^\n]", line) > 0) {
        if (line[0] == 't') {
            num_triangles++;
        }
        if (line[0] == 'c') {
            sscanf(line + 2, "%d%d%d", &r, &g, &b);
            kolorea = 1;
        }
    }
    fclose(obj_file);

    *hkopptr = num_triangles;
    *hptrptr = (hiruki *)malloc(num_triangles * sizeof(hiruki));
    if (*hptrptr == NULL) {
        printf("Error: no se pudo asignar memoria para los triángulos\n");
        return -1;
    }

    // Segunda lectura del archivo para llenar los datos de los triángulos
    obj_file = fopen(fitxiz, "r");
    i = 0;
    while (fscanf(obj_file, "\n%[^\n]", line) > 0) {
        if (line[0] == 't') {
            zkop = sscanf(line + 2, "%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f",
                          &v1, &v2, &v3, &v4, &v5, &v6, &v7, &v8, &v9,
                          &v10, &v11, &v12, &v13, &v14, &v15);
            printf("Número de valores en la línea = %d\n", zkop);
            if (zkop == 15) {
                (*hptrptr)[i].p1.x = v1;
                (*hptrptr)[i].p1.y = v2;
                (*hptrptr)[i].p1.z = v3;
                (*hptrptr)[i].p1.u = v4;
                (*hptrptr)[i].p1.v = v5;
                (*hptrptr)[i].p2.x = v6;
                (*hptrptr)[i].p2.y = v7;
                (*hptrptr)[i].p2.z = v8;
                (*hptrptr)[i].p2.u = v9;
                (*hptrptr)[i].p2.v = v10;
                (*hptrptr)[i].p3.x = v11;
                (*hptrptr)[i].p3.y = v12;
                (*hptrptr)[i].p3.z = v13;
                (*hptrptr)[i].p3.u = v14;
                (*hptrptr)[i].p3.v = v15;
            } else if (zkop == 9) {
                (*hptrptr)[i].p1.x = v1;
                (*hptrptr)[i].p1.y = v2;
                (*hptrptr)[i].p1.z = v3;
                (*hptrptr)[i].p2.x = v4;
                (*hptrptr)[i].p2.y = v5;
                (*hptrptr)[i].p2.z = v6;
                (*hptrptr)[i].p3.x = v7;
                (*hptrptr)[i].p3.y = v8;
                (*hptrptr)[i].p3.z = v9;
                // Set default texture coordinates
                (*hptrptr)[i].p1.u = (*hptrptr)[i].p1.v = 0.0;
                (*hptrptr)[i].p2.u = (*hptrptr)[i].p2.v = 0.0;
                (*hptrptr)[i].p3.u = (*hptrptr)[i].p3.v = 0.0;
            } else {
                printf("Error de formato en el archivo %s\n", fitxiz);
                free(*hptrptr);
                fclose(obj_file);
                return -1;
            }
            i++;
        }
    }
    fclose(obj_file);

    // Asignar el color si se especificó en el archivo
    if (kolorea) {
        *rgbptr = (unsigned char *)malloc(3 * sizeof(unsigned char));
        if (*rgbptr == NULL) {
            printf("Error: no se pudo asignar memoria para el color\n");
            free(*hptrptr);
            return -1;
        }
        (*rgbptr)[0] = (unsigned char)r;
        (*rgbptr)[1] = (unsigned char)g;
        (*rgbptr)[2] = (unsigned char)b;
        printf("koloriak kargatzen\n");

        return 9;  // Indicar que se cargaron los colores
    }

    return 15;  // Indicar que no se cargaron colores
}





/* * /
void main(int argc, char*argv[])
{
int num_triangles,i; 
hiruki *tptr;
unsigned char kolore[3];

printf("%s fitxategitik triangeluak kargatzera\n",argv[1]);
cargar_triangulos_color(argv[1], &num_triangles, &tptr,&(kolore[0]));
for (i=0; i<num_triangles; i++)
	{
	printf("t: (%.1f, %.1f, %.1f) (%.1f, %.1f)", tptr[i].p1.x, tptr[i].p1.y, tptr[i].p1.z, tptr[i].p1.u, tptr[i].p1.v);
	printf("   (%.1f, %.1f, %.1f) (%.1f, %.1f)", tptr[i].p2.x, tptr[i].p2.y, tptr[i].p2.z, tptr[i].p2.u, tptr[i].p2.v);
	printf("   (%.1f, %.1f, %.1f) (%.1f, %.1f)\n", tptr[i].p3.x, tptr[i].p3.y, tptr[i].p3.z, tptr[i].p3.u, tptr[i].p3.v);
	}
printf("kolorea: %d, %d,%d\n",kolore[0],kolore[1],kolore[2]);
}


 // */
