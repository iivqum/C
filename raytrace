#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <png.h>

int SCNW = 256;
int SCNH = 256;

typedef struct P3D{float x,y,z;}P3D;
typedef struct Ray{P3D p,d;}Ray;
typedef struct HitRecord{float t;P3D p,n;}HitRecord;
typedef struct Sphere{P3D p;float r;}Sphere;

void PSet(P3D *p,float x,float y,float z){p->x=x;p->y=y;p->z=z;}
void PSetP(P3D *p0,P3D *p1){p0->x=p1->x;p0->y=p1->y;p0->z=p1->z;}
void PAdd(P3D *p0,P3D *p1){p0->x+=p1->x;p0->y+=p1->y;p0->z+=p1->z;}
void PSub(P3D *p0,P3D *p1){p0->x-=p1->x;p0->y-=p1->y;p0->z-=p1->z;}
float PDot(P3D *p0,P3D *p1){
	return p0->x*p1->x+p0->y*p1->y+p0->z*p1->z;
}
void PScl(P3D *p0,float s){p0->x*=s;p0->y*=s;p0->z*=s;}

int RayHitSphere(Ray *r,Sphere *s)
{
	P3D oc;
	PSetP(&oc,&r->p);
	PSub(&oc,&s->p);
	float a = PDot(&r->d,&r->d);
	float b = 2*PDot(&oc,&r->d);
	float c = PDot(&oc,&oc)-(s->r*s->r);
	float d = b*b-4*a*c;
	if (d<=0)
		return 0;
	float d2 = sqrtf(d);
	float root = (-b-d2)/(2*a);
	if (root>0)return 1;
	root = (-b+d2)/(2*a);
	if (root>0)return 1;
	return 0;
}

struct Sphere *SphereInit(float x,float y,float z,float r)
{
	struct Sphere *sp = malloc(sizeof(struct Sphere));
	if (sp==NULL)return NULL;
	PSet(&sp->p,x,y,z);
	sp->r = r;
	return sp;
}

int main(int argc, char *argv[])
{
	FILE *fp = fopen(argv[1],"wb");
	if (fp==NULL)abort();
	png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,NULL,NULL,NULL);
	if (!png_ptr)abort();
	png_infop info_ptr = png_create_info_struct(png_ptr);
	if (!info_ptr){
		//png_destroy_write_struct(&png_ptr,NULL);
		abort();
	}
	if (setjmp(png_jmpbuf(png_ptr))){
		//png_destroy_write_struct(&png_ptr,&info_ptr);
		abort();
	}
	png_init_io(png_ptr,fp);
	png_set_IHDR(
		png_ptr,
		info_ptr,
		SCNW, SCNH,
		8,
		PNG_COLOR_TYPE_RGB,
		PNG_INTERLACE_NONE,
		PNG_COMPRESSION_TYPE_DEFAULT,
		PNG_FILTER_TYPE_DEFAULT
	);
	png_write_info(png_ptr,info_ptr);

	Sphere *s = SphereInit(0,0,-20,10);
	Ray r;
	PSet(&r.p,0,0,0);
	float ar = SCNW/SCNH;
	png_byte row[SCNW*3*sizeof(png_byte)];
	for (int i=0;i<SCNH;i++){
		for (int j=0;j<SCNW;j++){
			float u = (2*((float)j/((float)SCNW-1))-1)*ar;
			float v = -(2*((float)i/((float)SCNH-1))-1);
			PSet(&r.d,u,v,-1);
			if (RayHitSphere(&r,s)){
				row[j*3] = 255;
				row[j*3+1] = 255;
				row[j*3+2] = 255;
			}else{
				row[j*3] = 0;
				row[j*3+1] = 0;
				row[j*3+2] = 0;
			}
		}
		png_write_row(png_ptr,(png_bytep)row);
	}
	png_write_end(png_ptr,NULL);

	png_destroy_write_struct(&png_ptr,&info_ptr);
	fclose(fp);
	return 0;
}
