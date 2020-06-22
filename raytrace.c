#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <png.h>

int SCNW = 512;
int SCNH = 512;
int NSAMPLES = 1;

float EPSL = 1e-3;
float PI = 3.1415;

#define NSPHERES 2

typedef struct Color{unsigned char r,g,b;}Color;
typedef struct P3D{float x,y,z;}P3D;
typedef struct Ray{P3D p,d;}Ray;
typedef struct HitRecord{float t;P3D p,n;}HitRecord;
typedef struct Sphere{P3D p;float r;}Sphere;

Sphere spherelist[NSPHERES] = {
	{{0,0,-30},10},
	{{0,-1010,-30},1000},
};

void PSet(P3D *p,float x,float y,float z){p->x=x;p->y=y;p->z=z;}
void PSetP(P3D *p0,P3D *p1){p0->x=p1->x;p0->y=p1->y;p0->z=p1->z;}
void PAdd(P3D *p0,P3D *p1){p0->x+=p1->x;p0->y+=p1->y;p0->z+=p1->z;}
void PSub(P3D *p0,P3D *p1){p0->x-=p1->x;p0->y-=p1->y;p0->z-=p1->z;}
float PDot(P3D *p0,P3D *p1){
	return p0->x*p1->x+p0->y*p1->y+p0->z*p1->z;
}
void PScl(P3D *p,float s){p->x*=s;p->y*=s;p->z*=s;}
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

Sphere *SphereInit(float x,float y,float z,float r)
{
	struct Sphere *sp = malloc(sizeof(struct Sphere));
	if (sp==NULL)ErrNoMem();
	PSet(&sp->p,x,y,z);
	sp->r = r;
	return sp;
}

float RayHitSphere(Ray *r,Sphere *s)
{
	P3D oc;
	PSetP(&oc,&r->p);
	PSub(&oc,&s->p);
	float a = PDot(&r->d,&r->d);
	float b = 2*PDot(&oc,&r->d);
	//large radii might cause problems with float
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
	rec->t = FLT_MAX;
	Sphere *s = NULL;
	for (int i=0;i<NSPHERES;i++){
		float t2 = RayHitSphere(r,spherelist+i);
		if (t2<=0||t2>rec->t)
			continue;
		rec->t = t2;
		s = spherelist+i;
	}
	if (s==NULL)return 0;
	PSetP(&rec->p,&r->d);
	PScl(&rec->p,rec->t);
	PAdd(&rec->p,&r->p);
	//sphere normal (p-c)/||(p-c)||
	PSetP(&rec->n,&rec->p);
	PSub(&rec->n,&s->p);
	PScl(&rec->n,1/s->r);
	return 1;
}

P3D *TraceRay(Ray *r,int depth)
{
	if (depth<=0)return P3DInit(0,0,0);
	HitRecord rec;
	if (!GetNearestHit(r,&rec))
		return P3DInit(1,1,1);
	P3D d,s;
	RandomInSphere(&d);
	PAdd(&d,&rec.n);

	PSetP(&s,&rec.p);
	PAdd(&s,&d);
	PSub(&s,&rec.p);
	if (PDot(&s,&rec.n)<=0)
		return P3DInit(0,0,0);
	Ray r2;
	PSetP(&r2.p,&rec.p);
	PSetP(&r2.d,&s);
	P3D *col = TraceRay(&r2,depth-1);
	PScl(col,0.5);
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

	Ray r;
	PSet(&r.p,0,0,0);
	float ar = (float)SCNW/(float)SCNH;
	png_byte row[SCNW*3*sizeof(png_byte)];
	for (int i=0;i<SCNH;i++){
		for (int j=0;j<SCNW;j++){
			P3D avg;
			PSet(&avg,0,0,0);
			for (int k=0;k<NSAMPLES;k++){
				float u = (2*(((float)j+RandomFrac())/((float)(SCNW-1)))-1)*ar;
				float v = -(2*(((float)i+RandomFrac())/((float)(SCNH-1)))-1);
				PSet(&r.d,u,v,-1);
				P3D *col = TraceRay(&r,50);
				PAdd(&avg,col);
				free(col);
			}
			PScl(&avg,1/(float)NSAMPLES);
			avg.x = sqrtf(avg.x);
			avg.y = sqrtf(avg.y);
			avg.z = sqrtf(avg.z);
			if (avg.x>1)avg.x=1;
			if (avg.y>1)avg.y=1;
			if (avg.z>1)avg.z=1;
			row[j*3] = (png_byte)(avg.x*255);
			row[j*3+1] = (png_byte)(avg.y*255);
			row[j*3+2] = (png_byte)(avg.z*255);			
/*
			float u = (2*((float)j/((float)(SCNW-1)))-1)*ar;
			float v = -(2*((float)i/((float)(SCNH-1)))-1);
			PSet(&r.d,u,v,-1);
			P3D *col = TraceRay(&r,50);
			row[j*3] = (png_byte)(sqrtf(col->x)*255);
			row[j*3+1] = (png_byte)(sqrtf(col->y)*255);
			row[j*3+2] = (png_byte)(sqrtf(col->z)*255);
			free(col);
*/
		}
		png_write_row(png_ptr,(png_bytep)row);
	}
	png_write_end(png_ptr,NULL);
	png_destroy_write_struct(&png_ptr,&info_ptr);
	fclose(fp);
	return 0;
}
