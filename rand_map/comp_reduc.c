#include<stdio.h>
#define lh 6

main()
{
    FILE *f;
    int x,y,c,ni,no,ti,to,lkj=1,i,j;
    int w,h,g,*p;
    char ch,name[30];

    printf("input file : ");
    for(i=0;(name[i]=getchar())!='\n';i++); name[i]=0;
    char in[i+1];
    for(j=0;j<i+1;j++) in[j]=name[j];

    printf("\noutput file : ");
    for(i=0;(name[i]=getchar())!='\n';i++); name[i]=0;
    char out[i+1];
    for(j=0;j<i+1;j++) out[j]=name[j];

    printf("\n================================\n\n");

    f=fopen(in,"r");

    printf(in);
    printf("\n\n");

    for(c=0;c<lh;c++){
        while ((ch=getc(f))!='\n') printf("%c",ch);
        printf("\n");}

    fscanf(f,"%d %d %d",&w,&h,&g);
    printf("%d %d %d",w,h,g);

    ti=h*w;

    for(c=0,ni=0;c<ti;c++){
        fscanf(f,"%d",&x);
        if (x==0) ni++;}

    fclose(f);
    f=fopen(in,"r");

    for (c=0;c<lh+1;c++) while ((ch=getc(f))!='\n');

    for(c=0,y=0;c<h*w;c++){
        fscanf(f,"%d",&x);
        if (x>y) y=x;}

    fclose(f);

    printf("\n\nMax value = %d\n",y);
    printf("\nPixels = %d",ti);
    printf("\nZeros = %d",ni);
    printf("\nNonzero pixels = %d\n\n",ti-ni);

    f=fopen(out,"r");

    printf("================================\n\n");
    printf(out);
    printf("\n\n");

    for(c=0;c<lh;c++){
        while ((ch=getc(f))!='\n') printf("%c",ch);
        printf("\n");}

    fscanf(f,"%d %d %d",&w,&h,&g);
    printf("%d %d %d",w,h,g);

    to=h*w;

    for(c=0,no=0;c<to;c++){
        fscanf(f,"%d",&x);
        if (x==0) no++;}

    fclose(f);
    f=fopen(out,"r");

    for (c=0;c<lh+1;c++) while ((ch=getc(f))!='\n');

    for(c=0,y=0;c<h*w;c++){
        fscanf(f,"%d",&x);
        if (x>y) y=x;}

    fclose(f);

    printf("\n\nMax value = %d\n",y);
    printf("\nPixels = %d",to);
    printf("\nZeros = %d",no);
    printf("\nNonzero pixels = %d\n\n",to-no);
    printf("================================\n\n");


    printf("Calculation variation = %.1f%%\n\n",100*(float)((to-no)-(ti-ni))/(float)(ti-ni));
    printf("================================\n");
    getchar();

}



