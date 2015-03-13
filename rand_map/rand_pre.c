#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#define lh 6                /*Number of lines in the headers of the lumlp files before the line of the width, height and gain*/
#define seed 2              /*Seed of the random number generator*/
#define use_default_data 0  /*If not equal to zero, it will use default parameters without asking for it*/

void make_areas(char* nlumin)
{
    FILE *f=NULL;
    int ox,oy,x,y,z,w,h,g,i,j,it,jt,itt=0,jtt=0,pi,pj,t,n,k,lim,s,u;
    float r1,imp,d,p,mean;
    char ch,name[50];

    srand(seed);

/*Acquiring the values of the parameters for making the areas*/
    if (use_default_data==0)
    {
    printf("\nox = ");          /*x coordinate of the observation point starting with one and from the left*/
        scanf("%d",&ox);
    printf("\noy = ");          /*y coordinate of the observation point starting with one and from the bottom*/
        scanf("%d",&oy);
    printf("\nr1 = ");          /*Radius that determine where the probability function (witch is proportional to 1/r) is equal to 1*/
        scanf("%f",&r1);
    printf("\nimp = ");         /*Relative importance of the luminosity of each pixel to the probability of being choosed to represent an area*/
        scanf("%f",&imp);
    }

/*Default parameters. Will be active only if use_default_data is not defined as 0*/
    else
    {
        ox=302;
        oy=58;
        r1=1;
        imp=0;
    }

/*Open the file containing the luminosity of all the zones.*/
    sprintf(name,"%s_recombined.pgm",nlumin);
    f=fopen(name,"r");
/*Skip the headers*/
    while ((ch=getc(f))!='\n');
/*Read the width, height and gain*/
    fscanf(f,"%d %d %d",&w,&h,&g);
/*Create the required tables*/
    char m[h][w];   /*Will contain the values of the randmap. That is a table of 0 and 1 dermining witch pixel will be calculated in ILLUMINA.*/
    unsigned short int mi[h][w],mj[h][w];   /*Will contain the information of what pixel in m[][] each pixel will be merged. In other words, it is the information of the areas.*/
    unsigned short int m1[h][w];    /*Will contain the luminosity of each pixel*/
/*Transform the cartesian coordinates of the observation point in table coordinates*/
    pi=h-oy;
    pj=ox-1;
/*Put the values of the luminosity file in m1[][]*/
    mean=0.;
    for(i=0;i<h;i++)
        for(j=0;j<w;j++)
        {
            fscanf(f,"%d",&x);
            m1[i][j]=x;
            mean=mean+(float)m1[i][j];
        }
    fclose(f);
    mean=mean/((float)w*(float)h);
    printf("Mean value %.0f\n",mean);
/*Make the randmap*/
    for(i=0;i<h;i++)
        for(j=0;j<w;j++)
        {
            d=sqrt((pi-i)*(pi-i)+(pj-j)*(pj-j));    /*Distance to the observator*/
            p=(r1*r1)/(d*d)*imp*(float)m1[i][j]/mean;    /*Simplified 1/r2 function that pass through 1 at r1 and as a vertical asymptote at r=0*/
            if (((float)(rand()%32769)/32768)<p || p<0 || d<=r1 ) m[i][j]=1;    /*Put 1 in m[][] if a random float number between 0 and 1 is under the probability function.*/
            else m[i][j]=0;
        }
/*Save m[][] in a file*/
    f=fopen("randmap.pgm","w");
    fprintf(f,"P2\n%d %d %d\n",w,h,1);
    for(i=0;i<h;i++)
        for(j=0;j<w;j++)
            fprintf(f,"%d\n",m[i][j]);
    fclose(f);

/*Calculate the areas by determining for each pixel what is the nearest pixel equal to 1 in m[][]. Each pixel equal to one in m[][] will therefore be the center of an area containing all the pixels near to it.*/
    printf("\nGENERATING THE AREAS :\n");

    for(i=0;i<h;i++)
    {
        for(j=0;j<w;j++)
        {
            for(it=i,jt=j,n=1,z=(w+h)*(w+h),lim=w+h;n<=lim;n++)
                for(s=0;s<=1;s++)
                    for(k=0;k<n;k++)
                    {
                        if( it>=0 && it<h && jt>=0 && jt<w )
                            if(m[it][jt]==1)
                            {
                                y=(i-it)*(i-it)+(j-jt)*(j-jt);
                                if (y<z)
                                {
                                    if(lim==w+h)lim=2*sqrt(y)+1;
                                    z=y;
                                    itt=it;
                                    jtt=jt;
                                }
                            }
                        if(s==0)it+=pow(-1,n);
                        if(s==1)jt+=pow(-1,n);
                    }
            mi[i][j]=itt;   /*Put the information of the areas in mi and mj.*/
            mj[i][j]=jtt;
        }
        if (i%(h/100)==0) printf("%.0f%%\r",(float)i/h*100);    /*Print the percentage of the calculation that is done.*/
    }
    printf("%.0f%%\n\n",(float)i/h*100);

/*Save mi and mj in files.*/
    sprintf(name,"%s_mi.rand",nlumin);
    printf("Create %s\n",name);
    f=fopen(name,"w");
    for(i=0;i<h;i++)
        for(j=0;j<w;j++)
            fprintf(f,"%d\n",mi[i][j]);
    fclose(f);
    sprintf(name,"%s_mj.rand",nlumin);
    printf("Create %s\n",name);
    f=fopen(name,"w");
    for(i=0;i<h;i++)
        for(j=0;j<w;j++)
            fprintf(f,"%d\n",mj[i][j]);
    fclose(f);
}

void make_lumlp_out(char* nlumin)
{
    FILE *f=NULL,*fn=NULL,*fmi=NULL,*fmj=NULL;

    int mi,mj,w,h,g,x,i,j,nzone;
    char name[50],ch;

    for(nzone=1;nzone<=9;nzone++)
    {
        sprintf(name,"%s_0%d.pgm",nlumin,nzone);
        f=fopen(name,"r");

        if (f==0) break;

        sprintf(name,"%s_mi.rand",nlumin);
        fmi=fopen(name,"r");
        sprintf(name,"%s_mj.rand",nlumin);
        fmj=fopen(name,"r");
        sprintf(name,"%s_0%d_new.pgm",nlumin,nzone);
        fn=fopen(name,"w");

        printf("Create %s\n",name);

        for(x=0;x<lh;x++){
            while ((ch=getc(f))!='\n') fprintf(fn,"%c",ch);
                fprintf(fn,"\n");}

        fscanf(f,"%d %d %d",&w,&h,&g);

        int m[h][w];

        for(i=0;i<h;i++)
            for(j=0;j<w;j++)
                m[i][j]=0;

        for(i=0;i<h;i++)
            for(j=0;j<w;j++)
            {
                fscanf(fmi,"%d",&mi);
                fscanf(fmj,"%d",&mj);
                fscanf(f,"%d",&x);
                m[mi][mj]+=x;
            }

        fclose(fmi);
        fclose(fmj);
        fclose(f);

        for(i=g=0;i<h;i++)
            for(j=0;j<w;j++)
                if(m[i][j]>g) g=m[i][j];

        fprintf(fn,"%d %d %d\n",w,h,g);

        for(i=0;i<h;i++)
            for(j=0;j<w;j++)
                fprintf(fn,"%d\n",m[i][j]);

        fclose(fn);
    }

    sprintf(name,"%s_recombined.pgm",nlumin);
    f=fopen(name,"r");
    sprintf(name,"%s_mi.rand",nlumin);
    fmi=fopen(name,"r");
    sprintf(name,"%s_mj.rand",nlumin);
    fmj=fopen(name,"r");
    sprintf(name,"%s_recombined_new.pgm",nlumin);
    fn=fopen(name,"w");

    printf("Create %s\n",name);

    while ((ch=getc(f))!='\n') fprintf(fn,"%c",ch);
        fprintf(fn,"\n");

    fscanf(f,"%d %d %d",&w,&h,&g);

    int m[h][w];

    for(i=0;i<h;i++)
        for(j=0;j<w;j++)
            m[i][j]=0;

    for(i=0;i<h;i++)
        for(j=0;j<w;j++)
        {
            fscanf(fmi,"%d",&mi);
            fscanf(fmj,"%d",&mj);
            fscanf(f,"%d",&x);
            m[mi][mj]+=x;
        }

    fclose(fmi);
    fclose(fmj);
    fclose(f);

    for(i=g=0;i<h;i++)
        for(j=0;j<w;j++)
            if(m[i][j]>g) g=m[i][j];

    fprintf(fn,"%d %d %d\n",w,h,g);

    for(i=0;i<h;i++)
        for(j=0;j<w;j++)
            fprintf(fn,"%d\n",m[i][j]);

    fclose(fn);
}

void reflec(char* nlumin,char* nrefin)
{
    FILE *ffull=NULL,*fnew=NULL,*fp=NULL,*fref=NULL,*frefnew=NULL,*fmi=NULL,*fmj=NULL;
    int x,i,j,mi,mj,w,h,g,r;
    float p,pt;
    char name[50],ch;

    sprintf(name,"%s_mi.rand",nlumin);
    fmi=fopen(name,"r");
    sprintf(name,"%s_mj.rand",nlumin);
    fmj=fopen(name,"r");
    sprintf(name,"%s_recombined.pgm",nlumin);
    ffull=fopen(name,"r");
    sprintf(name,"%s_recombined_new.pgm",nlumin);
    fnew=fopen(name,"r");
    sprintf(name,"%s_p.rand",nlumin);
    printf("Create %s\n",name);
    fp=fopen(name,"w");
    printf("Create %s\n",name);
    sprintf(name,"%s_new.pgm",nrefin);
    frefnew=fopen(name,"w");
    printf("Create %s\n",name);

    while ((ch=getc(ffull))!='\n');
    fscanf(ffull,"%d %d %d",&w,&h,&g);

    while ((ch=getc(fnew))!='\n');
    fscanf(fnew,"%d %d %d",&w,&h,&g);

    int m[h][w];

    for(i=0;i<h;i++)
        for(j=0;j<w;j++)
            fscanf(fnew,"%d",&m[i][j]);

    fclose(fnew);

    for(i=0;i<h;i++)
        for(j=0;j<w;j++)
        {
            fscanf(fmi,"%d",&mi);
            fscanf(fmj,"%d",&mj);
            fscanf(ffull,"%d",&x);
            if (m[mi][mj]==0) p=0;
            else p=(float)x/m[mi][mj];
            fprintf(fp,"%f\n",p);
        }

    fclose(fmi);
    fclose(fmj);
    fclose(ffull);
    fclose(fp);

    sprintf(name,"%s.pgm",nrefin);
    printf("\nRead %s :\n\n",name);
    fref=fopen(name,"r");

    for(x=0;x<lh+1;x++){
        while ((ch=getc(fref))!='\n') {printf("%c",ch); fprintf(frefnew,"%c",ch);}
        printf("\n"); fprintf(frefnew,"%c",ch);}

    fscanf(fref,"%d %d %d",&w,&h,&g);
    printf("%d %d %d\n",w,h,g);
    fprintf(frefnew,"%d %d %d\n",w,h,g);

    float ref[h][w];

    sprintf(name,"%s_p.rand",nlumin);
    fp=fopen(name,"r");

    for(i=0;i<h;i++)
        for(j=0;j<w;j++)
            {
                fscanf(fp,"%f",&p);
                fscanf(fref,"%d",&r);
                ref[i][j]=p*r;
            }

    fclose(fref);
    fclose(fp);

    sprintf(name,"%s_mi.rand",nlumin);
    fmi=fopen(name,"r");
    sprintf(name,"%s_mj.rand",nlumin);
    fmj=fopen(name,"r");

    for(i=0;i<h;i++)
        for(j=0;j<w;j++)
            {
                fscanf(fmi,"%d",&mi);
                fscanf(fmj,"%d",&mj);
                ref[mi][mj]+=ref[i][j];
            }

    fclose(fmi);
    fclose(fmj);

    sprintf(name,"%s_mi.rand",nlumin);
    fmi=fopen(name,"r");
    sprintf(name,"%s_mj.rand",nlumin);
    fmj=fopen(name,"r");

    for(i=0;i<h;i++)
        for(j=0;j<w;j++)
            {
                fscanf(fmi,"%d",&mi);
                fscanf(fmj,"%d",&mj);
                ref[i][j]=ref[mi][mj];
            }

    fclose(fmi);
    fclose(fmj);

    for(i=0,pt=0;i<h;i++)
        for(j=0;j<w;j++){
            p=ref[i][j];
            if (p>pt) pt=p;}
    g=(int)(pt+1);

    for(i=0;i<h;i++)
        for(j=0;j<w;j++)
            fprintf(frefnew,"%d\n",(int)(ref[i][j]));

    fclose(frefnew);
}

void recombine(char* nlumin)
{
    FILE *f=NULL;

    int x,i,j,nz,v,w,h,g;
    char nend[50],name[50],ch;

    sprintf(nend,"_0%d.pgm",1);

    sprintf(name,"%s%s",nlumin,nend);

    f=fopen(name,"r");

    if (f==0) {printf("\nERROR : %s not found",name); getchar();}

    for(x=0;x<lh;x++)while((ch=getc(f))!='\n');

    fscanf(f,"%d %d %d",&w,&h,&g);

    fclose(f);

    int m[h][w];

    for(i=0;i<h;i++)
        for(j=0;j<w;j++)
            m[i][j]=0;

    for(nz=1;nz<=9;nz++)
    {
        sprintf(nend,"_0%d.pgm",nz);

        sprintf(name,"%s%s",nlumin,nend);

        f=fopen(name,"r");

        if (f==0) break;

        printf("\nRead %s :\n\n",name);

        for(x=0;x<lh;x++){
            while ((ch=getc(f))!='\n') printf("%c",ch);
                printf("\n");}

        fscanf(f,"%d %d %d",&w,&h,&g);
        printf("%d %d %d\n",w,h,g);

        for(i=0;i<h;i++)
            for(j=0;j<w;j++)
            {
                fscanf(f,"%d",&v);
                m[i][j]+=v;
            }

        fclose(f);
    }

    sprintf(name,"%s_recombined.pgm",nlumin);
    printf("\nCreate %s\n",name);
    f=fopen(name,"w");

    fprintf(f,"P2\n%d %d %d\n",w,h,g);

    for(i=0;i<h;i++)
        for(j=0;j<w;j++)
            fprintf(f,"%d\n",m[i][j]);

    fclose(f);
}

int main()
{
    int i;
    char nlumin[50],nrefin[50];

/*Acquiring the name of the luminosity and reflectance map that will be modified by this program*/
    printf("lumlp input : ");
        for(i=0;(nlumin[i]=getchar())!='\n';i++); nlumin[i]='\0';
    printf("reflect input : ");
        for(i=0;(nrefin[i]=getchar())!='\n';i++); nrefin[i]='\0';

/*Recombine all the zones into one file*/
    recombine(nlumin);

/*Determine the areas*/
    make_areas(nlumin);

/*Make the new luminosity maps*/
    make_lumlp_out(nlumin);

/*Make the new reflectance map*/
    reflec(nlumin,nrefin);
    return 0;
}
