#include "gmp.h"
#include <stdio.h>
		#include <time.h>
void main()
{
    mpz_t a, b, c, r, a_copy;
    mpz_init(a);
    mpz_init(b);
    mpz_init(c);
    mpz_init(r);
    mpz_set_ui(a,5);
    mpz_set_ui(b,2);
    mpz_fdiv_qr(c,r,a,b);
    gmp_printf("c=a/b=%Zd/%Zd=%Zd...%Zd      ",a,b,c,r);
    mpz_cdiv_qr(c,r,a,b);
    gmp_printf("c=a/b=%Zd/%Zd=%Zd...%Zd\n",a,b,c,r);
    mpz_set_si(a,-5);
    mpz_set_ui(b,2);
    mpz_fdiv_qr(c,r,a,b);
    gmp_printf("c=a/b=%Zd/%Zd=%Zd...%Zd      ",a,b,c,r);
    mpz_cdiv_qr(c,r,a,b);
    gmp_printf("c=a/b=%Zd/%Zd=%Zd...%Zd\n",a,b,c,r);
    mpz_set_ui(a,5);
    mpz_set_si(b,-2);
    mpz_fdiv_qr(c,r,a,b);
    gmp_printf("c=a/b=%Zd/%Zd=%Zd...%Zd      ",a,b,c,r);
    mpz_cdiv_qr(c,r,a,b);
    gmp_printf("c=a/b=%Zd/%Zd=%Zd...%Zd\n",a,b,c,r);

    int t[10]= {1,2,3,3,5,7,8,9,3,7};
    int tt[10]={2,4,1,3,1,6,3,6,5,8};
    int max = 10;
    int i;
    for(i = 0;i<max;i++)
    {
        if (tt[i]>t[i])
	{
		printf("yeah i=%d, max = %d\n",i,max);
		max--;
	}
    }
	printf("yeah i=%d, max = %d\n",i,max);

	clock_t start, end;
	double t_fdiv=0, t_cdiv=0, t_tdiv=0, t_divex;

	mpz_set_ui(a,23516871546525);
	mpz_mul(a,a,a);
	mpz_mul(a,a,a);
	mpz_set_si(b,2164466);
	mpz_mul(b,b,b);
	mpz_mul(b,b,b);
	mpz_mul(a,a,b);
	for (int iter=0;iter<100;iter++)
	{
		start = clock();
		for(i=0;i<100000;i++)
		{
			mpz_cdiv_q(c,a,b);
		}
		end = clock();
		t_cdiv += ((double) (end - start)) / CLOCKS_PER_SEC;
		start = clock();
		for(i=0;i<100000;i++)
		{
			mpz_tdiv_q(c,a,b);
		}
		end = clock();
		t_tdiv += ((double) (end - start)) / CLOCKS_PER_SEC;
		start = clock();
		for(i=0;i<100000;i++)
		{
			mpz_fdiv_q(c,a,b);
		}
		end = clock();
		t_fdiv += ((double) (end - start)) / CLOCKS_PER_SEC;
		start = clock();
		for(i=0;i<100000;i++)
		{
			mpz_divexact(c,a,b);
		}
		end = clock();
		t_divex += ((double) (end - start)) / CLOCKS_PER_SEC;
	}
	printf("time for fdiv: %f, cdiv:%f, tdiv:%f, divex:%f\n",t_fdiv, t_cdiv, t_tdiv, t_divex);
}
