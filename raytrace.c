#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <png.h>

int SCNW = 512;
int SCNH = 512;

float EPSL = 1e-3;
float PI = 3.1415;

typedef struct Color{unsigned char r,g,b;}Color;
typedef struct P3D{float x,y,z;}P3D;
typedef struct Ray{P3D p,d;}Ray;
typedef struct HitRecord{float t;P3D p,n;int hit;}HitRecord;
typedef struct Sphere{P3D p;float r;}Sphere;

void PSet(P3D *p,float x,float y,float z){p->x=x;p->y=y;p->z=z;}
void PSetP(P3D *p0,P3D *p1){p0->x=p1->x;p0->y=p1->y;p0->z=p1->z;}
void PAdd(P3D *p0,P3D *p1){p0->x+=p1->x;p0->y+=p1->y;p0->z+=p1->z;}
void PSub(P3D *p0,P3D *p1){p0->x-=p1->x;p0->y-=p1->y;p0->z-=p1->z;}
float PDot(P3D *p0,P3D *p1){
	return p0->x*p1->x+p0->y*p1->y+p0->z*p1->z;
}
void PScl(P3D *p0,float s){p0->x*=s;p0->y*=s;p0->z*=s;}

void ErrNoMem(void){perror("NOMEM");abort();}

float RandomFrac(void){
	return (float)rand()/(float)(RAND_MAX+1);
}

float RandomFracSign(void){
	return 2*RandomFrac()-1;
}

void RandomInSphere(P3D *p)
{
	float t = RandomFrac()*2*PI;
	float z = RandomFracSign();
	float r = sqrtf(1-z*z);
	PSet(p,r*cosf(t),r*sinf(t),z);
}

P3D *P3DInit(float x,float y,float z){
	P3D *ptr = malloc(sizeof(P3D));
	if (ptr==NULL)ErrNoMem();
	PSet(ptr,x,y,z);
	return ptr;
}

Color *ColorInit(unsigned char r,unsigned char g,unsigned char b){
	Color *ptr = malloc(sizeof(Color));
	if (ptr==NULL)ErrNoMem();
	ptr->r = r;
	ptr->g = g;
	ptr->b = b;
	return ptr;
}

Sphere *SphereInit(float x,float y,float z,float r)
{
	struct Sphere *sp = malloc(sizeof(struct Sphere));
	if (sp==NULL)ErrNoMem();
	PSet(&sp->p,x,y,z);
	sp->r = r;
	return sp;
}

#define NSPHERES 2

Sphere spherelist[NSPHERES];

float RayHitSphere(Ray *r,Sphere *s)
{
	P3D oc;
	PSetP(&oc,&r->p);
	PSub(&oc,&s->p);
	float a = PDot(&r->d,&r->d);
	float b = 2*PDot(&oc,&r->d);
	float c = PDot(&oc,&oc)-(s->r*s->r);
	float d = b*b-4*a*c;
	if (d<=0)return 0;
	float d2 = sqrtf(d);
	float root = (-b-d2)/(2*a);
	if (root>EPSL)return root;
	root = (-b+d2)/(2*a);
	if (root>EPSL)return root;
	return -1;
}

int GetNearestHit(Ray *r,HitRecord *rec)
{
	float t = FLT_MAX;
	Sphere *s = NULL;
	for (int i=0;i<NSPHERES;i++){
		float t2 = RayHitSphere(r,spherelist+i);
		if (t2>0&&t2<t){
			t = t2;
			s = spherelist+i;
		}
	}
	if (s==NULL)return 0;
	PSetP(&rec->p,&r->d);
	PScl(&rec->p,rec->t);
	PAdd(&rec->p,&r->p);
	//sphere normal (p-c)/||(p-c)||
	PSetP(&rec->n,&rec->p);
	PSub(&rec->n,&s->p);
	PScl(&rec->n,PDot(&rec->n,&rec->n));
	return 1;
}

Color *TraceRay(Ray *r,int depth)
{
	if (depth<=0)return ColorInit(0,0,0);
	HitRecord rec;
	if (!GetNearestHit(r,&rec))
		return ColorInit(255,255,255);
	P3D d,s;
	RandomInSphere(&d);
	PAdd(&d,&rec.n);

	PSetP(&s,&rec.p);
	PAdd(&s,&d);

	PSub(&s,&rec.p);
	Ray r2;
	PSetP(&r2.p,&rec.p);
	PSetP(&r2.d,&s);
	Color *col = TraceRay(&r2,depth-1);
	col->r = (float)col->r*(float)0.5;
	col->g = (float)col->g*0.5;
	col->b = (float)col->b*0.5;
	printf("%d,%d,%d\n",col->r,col->g,col->b);
	return col;
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

	PSet(&spherelist[0].p,0,0,-30);
	spherelist[0].r = 10;

	PSet(&spherelist[1].p,0,-1010,-30);
	spherelist[1].r = 1000;

	Ray r;
	PSet(&r.p,0,0,0);
	float ar = SCNW/SCNH;
	png_byte row[SCNW*3*sizeof(png_byte)];
	for (int i=0;i<SCNH;i++){
		for (int j=0;j<SCNW;j++){
			float u = (2*((float)j/((float)SCNW-1))-1)*ar;
			float v = -(2*((float)i/((float)SCNH-1))-1);
			PSet(&r.d,u,v,-1);
			Color *col = TraceRay(&r,50);
			//printf("%f,%f,%f\n",col->x,col->y,col->z);
			row[j*3] = (png_byte)col->r;
			row[j*3+1] = (png_byte)col->g;
			row[j*3+2] = (png_byte)col->b;
			free(col);
		}
		png_write_row(png_ptr,(png_bytep)row);
	}
	png_write_end(png_ptr,NULL);

	png_destroy_write_struct(&png_ptr,&info_ptr);
	fclose(fp);
	return 0;
}
