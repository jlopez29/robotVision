#include <stdio.h>
#include <math.h>
#define PICSIZE 256
#define MAXMASK 100

         int    pic[PICSIZE][PICSIZE];
         double outpicx[PICSIZE][PICSIZE];
         double outpicy[PICSIZE][PICSIZE];
         int    edgeflag[PICSIZE][PICSIZE];
         int    peaks[PICSIZE][PICSIZE];
         int    final[PICSIZE][PICSIZE];
         double maskx[MAXMASK][MAXMASK];
         double masky[MAXMASK][MAXMASK];
         double convx[PICSIZE][PICSIZE];
         double convy[PICSIZE][PICSIZE];
         double ival[256][256], mag[256][256];


void recurse(int i, int j,int lo);

main()
{

    FILE *fo1, *fo2, *fo3, *fi;
    int     i,j,p,q,s,t,mr,centx,centy;
    double  maskvalx,maskvaly,sumx,sumy,sig=1.0,maxival,minival,maxval,perc;

    fi=fopen("garb34","rb");

    fo1=fopen("cannygradient.pgm","wb");
    fo2=fopen("cannypeaks.pgm","wb");
    fo3=fopen("cannyfinal.pgm","wb");

    /*  Creating Headers for output files.  */
    fprintf(fo1,"P5\n");
    fprintf(fo1,"%d %d\n", 256, 256);
    fprintf(fo1,"255\n");

    fprintf(fo2,"P5\n");
    fprintf(fo2,"%d %d\n", 256, 256);
    fprintf(fo2,"255\n");

    fprintf(fo3,"P5\n");
    fprintf(fo3,"%d %d\n", 256, 256);
    fprintf(fo3,"255\n");


    /*   Copy input picture to pic array  */
    for (i=0;i<256;i++)
        { for (j=0;j<256;j++)
                {
                  pic[i][j]  =  getc (fi);
                }
        }


    printf("Enter a value for sigma.\n");
    scanf("%lf",&sig);

    printf("Enter histogram percentage:\n");
    scanf("%lf",&perc);


    mr = (int)(sig * 3);
    centx = (MAXMASK / 2);
    centy = (MAXMASK / 2);


    /*  Generate x and y masks  */
    for (p=-mr;p<=mr;p++)
        {  for (q=-mr;q<=mr;q++)
           {
               maskvalx = (-q/(sig*sig))*((exp(-1*(((p*p)+(q*q))/(2*(sig*sig))))));
               maskvaly = (-p/(sig*sig))*((exp(-1*(((p*p)+(q*q))/(2*(sig*sig))))));
              (maskx[p+centy][q+centx]) = maskvalx;
              (masky[p+centy][q+centx]) = maskvaly;
           }
        }

    /*  Get 2 sums for convolution */

    for (i=mr;i<=255-mr;i++)
        { for (j=mr;j<=255-mr;j++)
          {
             sumx = 0;
             sumy = 0;
             for (p=-mr;p<=mr;p++)
             {
                for (q=-mr;q<=mr;q++)
                {
                   sumx += pic[i+p][j+q] * maskx[p+centy][q+centx];
                   sumy += pic[i+p][j+q] * masky[p+centy][q+centx];
                }
             }
             outpicx[i][j] = sumx;
             outpicy[i][j] = sumy;
             convx[i][j] = sumx;
             convy[i][j] = sumy;
          }
        }

        /* Compute magnitude and check maximum (from Sobel) */

        maxival = 0;
        for (i=mr;i<256-mr;i++)
        { for (j=mr;j<256-mr;j++)
          {
             ival[i][j]=sqrt((double)((outpicx[i][j]*outpicx[i][j]) +
                                      (outpicy[i][j]*outpicy[i][j])));
             if (ival[i][j] > maxival)
                maxival = ival[i][j];

           }
        }

        double slope;

        for(i=mr;i<256-mr;i++){
           for(j=mr;j<256-mr;j++){

              if((convx[i][j]) == 0.0) {
                 convx[i][j] = .00001;
               }
              slope = convy[i][j]/convx[i][j];

              if( (slope <= .4142)&&(slope > -.4142)){
                 if((ival[i][j] > ival[i][j-1])&&(ival[i][j] > ival[i][j+1])){
                    peaks[i][j] = 1;
                  }
               }
              else if( (slope <= 2.4142)&&(slope > .4142)){
                 if((ival[i][j] > ival[i-1][j-1])&&(ival[i][j] > ival[i+1][j+1])){
                     peaks[i][j] = 1;
                   }
               }
              else if( (slope <= -.4142)&&(slope > -2.4142)){
                 if((ival[i][j] > ival[i+1][j-1])&&(ival[i][j] > ival[i-1][j+1])){
                     peaks[i][j] = 1;
                   }
               }else{
                 if((ival[i][j] > ival[i-1][j])&&(ival[i][j] > ival[i+1][j])){
                     peaks[i][j] = 1;
                  }
                }
             }
          }


        for(i=0;i<256;i++)
        { for(j=0;j<256;j++)
            {
                mag[i][j]=ival[i][j];
            }
        }

         /* Scale to 0~255 and print to outputgradient.pgm */
        for (i=0;i<256;i++)
          { for (j=0;j<256;j++)
            {
             ival[i][j] = (ival[i][j] / maxival) * 255;
             fprintf(fo1,"%c",(char)((int)(ival[i][j])));

            }
          }
        fclose(fo1);

        int x2= 255;
        int x=0;

        /* Print out peaks */

        for (i=0;i<256;i++)
          { for (j=0;j<256;j++)
            {
             if(peaks[i][j] == 1)
             fprintf(fo2,"%c",(char)((int)(x2)));
             else
             fprintf(fo2,"%c",(char)((int)(x)));

            }
          }

        int hi;
        int lo;
        int hist[256];
        /* Initialize Histogram */
        for(i=0; i < 256 ; i++)
        {
            hist[i]=0;
        }

        /* Populate Histogram */
        for(i=0; i <256 ; i++)
        {
            for (j=0; j < 256;j++)
            {
               hist[(int)ival[i][j]]++;
            }
        }


        double area;
        int sum=0;
        /* Calculate area */
        area = (perc/100) * (256*256);
        /* Set hi/lo */
        for(t=255; t > 0 ; t--)
        {
            sum+=hist[t];
            if(sum > (int)area)
               break;
        }

        hi=t;
        lo = hi*0.35;

        printf("Hi = %d\n",hi);
        printf("Lo = %d\n",lo);
        system("PAUSE");

        /* Chaining Peaks */
        for (i=0; i<256;i++)
        {
            for(j=0;j<256;j++)
            {
                if(peaks[i][j] == 1)
                {
                  if(ival[i][j] > hi)
                  {
                      peaks[i][j] = 0;
                      final[i][j] = 1;
                      recurse(i,j,lo);
                  }
                  else if(ival[i][j] < lo)
                  {
                      peaks[i][j] = 0;
                      final[i][j] = 0;
                  }
                }
            }
        }


    /* Print chaining peaks */
    for (i=0;i<256;i++)
    {
        for (j=0;j<256;j++)
        {
             if(final[i][j] == 1)
                fprintf(fo3,"%c",(char)((int)(x2)));
             else
                fprintf(fo3,"%c",(char)((int)(x)));

         }
     }
}

void recurse(int i, int j, int lo)
{

    int p,q;

    for(p=-1; p <= 1 ; p++)
    { for(q=-1; q <= 1; q++)
        {
            if(peaks[i+p][j+q]==1)
            {
                peaks[i+p][j+q]=0;
                final[i+p][j+q]=1;
                recurse(i+p,j+q,lo);

            }
        }
    }

}



