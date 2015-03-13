#include<stdio.h>
#define lh 1

main()
{
    FILE *mi,*mj,*f;
    int w=395,h=290,g,i,j,vi,vj;
    int rli[h],rlj[w],a[h][w];
    char nlumin[50],name[50];
    g=w*h;

    printf("lumlp input : ");
        for(i=0;(nlumin[i]=getchar())!='\n';i++); nlumin[i]=0;

    for(i=0;i<h;i++)
        rli[i]=rand()%(h+1);
    for(j=0;j<w;j++)
        rlj[j]=rand()%(w+1);

    sprintf(name,"%s_mi.rand",nlumin);
    mi=fopen(name,"r");
    sprintf(name,"%s_mj.rand",nlumin);
    mj=fopen(name,"r");

    for(i=0;i<h;i++)
        for(j=0;j<w;j++)
        {
            fscanf(mi,"%d",&vi);
            fscanf(mj,"%d",&vj);
            a[i][j]=rli[vi]*rlj[vj];
        }

    fclose(mi);
    fclose(mj);

    f=fopen("areas_rand.pgm","w");
    fprintf(f,"P2\n%d %d %d\n",w,h,g);

    for(i=0;i<h;i++)
        for(j=0;j<w;j++)
            fprintf(f,"%d\n",a[i][j]);

    fclose(f);

}

