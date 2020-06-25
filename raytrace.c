#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <png.h>

int SCNW = 256;
int SCNH = 256;
int NSAMPLES = 100;

float EPSL = 1e-3;
float PI = 3.1415;

typedef enum {DIFFUSE,METALLIC}material_type;

typedef struct Color{unsigned char r,g,b;}Color;
typedef struct P3D{float x,y,z;}P3D;
typedef struct Ray{P3D p,d;}Ray;
typedef struct Sphere{P3D p,color;material_type mt;float r;}Sphere;
typedef struct Plane{P3D p,color;P3D n;}Plane;
typedef struct HitRecord{float t;P3D p,n,*color;material_type mt;}HitRecord;

Sphere spherelist[] = {
	{{0,0,-30},{0.85,0.85,0.85},METALLIC,10},
	{{-20,0,-30},{0.85,0.85,0.85},DIFFUSE,10},
	{{20,0,-30},{0.85,0.85,0.85},DIFFUSE,10},
	{{0,-1010,-30},{0.85,0.85,0.85},DIFFUSE,1000}
};

Plane planelist[] = {

};

#define NSPHERES sizeof(spherelist)/sizeof(Sphere)
#define NPLANES sizeof(planelist)/sizeof(Sphere)

void PSet(P3D *p,float x,float y,float z){p->x=x;p->y=y;p->z=z;}
void PSetP(P3D *p0,P3D *p1){p0->x=p1->x;p0->y=p1->y;p0->z=p1->z;}
void PAdd(P3D *p0,P3D *p1){p0->x+=p1->x;p0->y+=p1->y;p0->z+=p1->z;}
void PSub(P3D *p0,P3D *p1){p0->x-=p1->x;p0->y-=p1->y;p0->z-=p1->z;}
void PMul(P3D *p0,P3D *p1){p0->x*=p1->x;p0->y*=p1->y;p0->z*=p1->z;}
float PDot(P3D *p0,P3D *p1){
	return p0->x*p1->x+p0->y*p1->y+p0->z*p1->z;
}
void PScl(P3D *p,float s){p->x*=s;p->y*=s;p->z*=s;}
void ErrNoMem(void){perror("NOMEM");abort();}

float RandomFrac(void){
	return (float)rand()/(float)(RAND_MAX+1);
}

void RandomInSphere(P3D *p)
{
	float t = RandomFrac()*2*PI;
	float z = 2*RandomFrac()-1;;
	float r = sqrtf(1-z*z);
	PSet(p,r*cosf(t),r*sinf(t),z);
}

P3D *P3DInit(float x,float y,float z){
	P3D *ptr = malloc(sizeof(P3D));
	if (ptr==NULL)ErrNoMem();
	PSet(ptr,x,y,z);
	return ptr;
}

int RayHitSphere(Ray *r,Sphere *s,float tmax,float *tout)
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
	if (root>EPSL&&root<tmax){*tout = root;return 1;}
	root = (-b+d2)/(2*a);
	if (root>EPSL&&root<tmax){*tout = root;return 1;}
	return 0;
}

int GetNearestHit(Ray *r,HitRecord *rec)
{
	rec->t = FLT_MAX;
	for (int i=0;i<NSPHERES;i++){
		Sphere *s = spherelist+i;
		if (!RayHitSphere(r,s,rec->t,&rec->t))
			continue;
		PSetP(&rec->p,&r->d);
		PScl(&rec->p,rec->t);
		PAdd(&rec->p,&r->p);
		//sphere normal
		PSetP(&rec->n,&rec->p);
		PSub(&rec->n,&s->p);
		PScl(&rec->n,1/s->r);	

		rec->color = &s->color;	
		rec->mt = s->mt;	
	}
	return rec->t<FLT_MAX;
}

int ComputeDiffuseRay(Ray *r,HitRecord *rec,Ray *out)
{
	P3D d,s;
	RandomInSphere(&d);
	PAdd(&d,&rec->n);
	PSetP(&s,&rec->p);
	PAdd(&s,&d);
	PSub(&s,&rec->p);
	if (PDot(&s,&rec->n)<=0)
		return 0;
	PSetP(&out->p,&rec->p);
	PSetP(&out->d,&s);
	return 1;
}

int ComputeReflectionRay(Ray *r,HitRecord *rec,Ray *out)
{
	P3D n;
	PSetP(&n,&rec->n);
	PSetP(&out->p,&rec->p);
	PSetP(&out->d,&r->d);	
	PScl(&n,2*PDot(&r->d,&rec->n));
	PSub(&out->d,&n);
	return 1;
}

P3D *TraceRay(Ray *r,int depth)
{
	if (depth<=0)return P3DInit(0,0,0);
	HitRecord rec;
	if (!GetNearestHit(r,&rec))
		return P3DInit(1,1,1);
	Ray new;
	switch(rec.mt){
		case METALLIC:
			ComputeReflectionRay(r,&rec,&new);
			break;
		case DIFFUSE:
			if (!ComputeDiffuseRay(r,&rec,&new))
				return P3DInit(0,0,0);
			break;
		//default:
		//	return P3DInit(0,0,0);
	}
	P3D *col = TraceRay(&new,depth-1);
	PMul(col,rec.color);
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
	//png_set_compression_level(png_ptr,0);
	png_write_info(png_ptr,info_ptr);
	Ray r;
	P3D avg;
	PSet(&r.p,0,0,0);
	float ar = (float)SCNW/(float)SCNH;
	png_byte row[SCNW*3*sizeof(png_byte)];
	for (int i=0;i<SCNH;i++){
		for (int j=0;j<SCNW;j++){
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
		printf("Progress:%d%c\r",(int)((float)i/(float)(SCNH-1)*100),'%');
		png_write_row(png_ptr,(png_bytep)row);
	}
	png_write_end(png_ptr,NULL);
	png_destroy_write_struct(&png_ptr,&info_ptr);
	fclose(fp);
	return 0;
}
