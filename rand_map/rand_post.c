#include<stdio.h>
#define lh 6

main()
{
    FILE *data,*mi,*mj,*f,*fw;
    int x,y,w,h,g,i,j,im,jm;
    float p;
    char ch,nlumpl[150],in[150],name[150];

    printf(".rand root file name (without extension) : ");
        for(i=0;(nlumpl[i]=getchar())!='\n';i++); nlumpl[i]=0;

    printf("\ninput root file (without extension) : ");
        for(i=0;(in[i]=getchar())!='\n';i++); in[i]=0;

    sprintf(name,"%s_p.rand",nlumpl);

    printf("\n%s\n",name);

    data=fopen(name,"r");

    sprintf(name,"%s_mi.rand",nlumpl);

    printf("%s\n",name);

    mi=fopen(name,"r");

    sprintf(name,"%s_mj.rand",nlumpl);

    printf("%s\n",name);

    mj=fopen(name,"r");

    sprintf(name,"%s.pgm",in);
    f=fopen(name,"r");
    printf("\n%s",in);
    printf("\n\n");

    sprintf(name,"%s_new.pgm",in);
    fw=fopen(name,"w");

    for(x=0;x<lh;x++){
        while ((ch=getc(f))!='\n') {printf("%c",ch);fprintf(fw,"%c",ch);}
        {printf("\n");fprintf(fw,"\n");}}

    fscanf(f,"%d %d %d\n",&w,&h,&g);
    printf("width = %d\nheight = %d\nmax = %d\n",w,h,g);

    unsigned int m1[h][w],m2[h][w];

    for(i=0;i<h;i++)
        for(j=0;j<w;j++)
            m1[i][j]=0;

    for(i=0,x=0;i<h;i++)
        for(j=0;j<w;j++)
        {
            fscanf(f,"%d",&x);
            fscanf(mi,"%d",&im);
            fscanf(mj,"%d",&jm);
            m1[im][jm]+=x;
        }

    fclose(f);
    fclose(mi);
    fclose(mj);

    sprintf(name,"%s_mi.rand",nlumpl);
    mi=fopen(name,"r");
    sprintf(name,"%s_mj.rand",nlumpl);
    mj=fopen(name,"r");

    for(i=0;i<h;i++)
        for(j=0;j<w;j++)
        {
            fscanf(mi,"%d",&im);
            fscanf(mj,"%d",&jm);
            fscanf(data,"%f",&p);
            m2[i][j]=(float)p*(float)m1[im][jm]+0.5;
        }

    fclose(data);
    fclose(mi);
    fclose(mj);

    for(i=g=0;i<h;i++)
        for(j=0;j<w;j++)
            if (m2[i][j]>g) g=m2[i][j];

    fprintf(fw,"%d %d %d\n",w,h,g);

    for(i=0;i<h;i++)
        for(j=0;j<w;j++)
            fprintf(fw,"%d\n",m2[i][j]);

    fclose(fw);
}
