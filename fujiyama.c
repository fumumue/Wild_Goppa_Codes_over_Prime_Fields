#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <math.h>
#include <assert.h>
#include <x86intrin.h> // SIMD命令を使用するためのヘッダファイル

#include "global-p.h"
#include "struct.h"
#include "chash-p.c"

#define SEPARABLE 0
#define MATRIX_SIZE 16
#define SHM_KEY 128

unsigned short g[K + 1] = {0};

// ランダム多項式の生成
static void
ginit(void)
{
    int j, count = 0, k = 0;
    unsigned short gg[K + 1] = {0};

    printf("in ginit\n");

    g[K] = 1;          // xor128();
    g[0] = rand() % N; // or N
    k = rand() % (K - 1);
    if (k > 0)
    {
        while (count < k)
        {
            printf("in whule\n");
            j = rand() % (K);
            if (j < K && j > 0 && g[j] == 0)
            {
                g[j] = rand() % N; // or N;
                count++;
            }
        }
    }

    for (j = 0; j < K + 1; j++)
        gg[j] = g[K - j];

    memcpy(g, gg, sizeof(g));
}

// OP型からベクトル型への変換
vec o2v(OP f)
{
    vec a = {0};
    int i;

    for (i = 0; i < K * E; i++)
    {
        if (f.t[i].a > 0 && f.t[i].n < K * E)
            a.x[f.t[i].n] = f.t[i].a;
    }

    return a;
}

// ベクトル型からOP型への変換
OP v2o(vec a)
{
    int i, j = 0;
    OP f = {0};

    // #pragma omp parallel for
    for (i = 0; i < K * E; i++)
    {
        if (a.x[i] > 0)
        {
            f.t[j].n = i;
            f.t[j++].a = a.x[i];
        }
    }

    return f;
}

unsigned short oinv(short a, unsigned short n)
{
    unsigned short i;

    if (a == 0)
        return 0;
    if (a < 0)
    {
        //printf("a=%d", a);
        a = N + a%N;
        //printf("-a=%d\n", a);
        // exit(1);
    }
    // if (a == 1)
    //     return 1;
    for (i = 1; i < n; i++)
    {
        if ((i * a) % N == 1)
            return i;
    }
    printf("no return\n");
    exit(1);
}

// 多項式の次数(default)
int deg(vec a)
{
    int i, n = 0, flg = 0;

    // #pragma omp parallel for
    for (i = 0; i < DEG; i++)
    {
        if (a.x[i] > 0)
        {
            n = i;
            flg = 1;
        }
    }
    if (flg == 0)
        return 0;

    return n;
}

// 多項式を表示する(default)
void printpol(vec a)
{
    int i, n;

    n = deg(a);

    // printf ("baka\n");
    //  assert(("baka\n", n >= 0));

    for (i = n; i > -1; i--)
    {
        if (a.x[i] > 0)
        {
            printf("%u*", a.x[i]);
            // if (i > 0)
            printf("x^%d", i);
            // if (i > 0)
            printf("+");
        }
    }
    //  printf("\n");

    return;
}

vec kof2(unsigned short c, vec f)
{
    int i, k;
    vec b = {0}, h = {0};

    c = oinv(c, N);
    printf("c=%d\n", c);
    // exit(1);
    b = f; // o2v(f);
    k = deg(b);
    printpol(b);
    printf(" =b debugi\n");
    for (i = 0; i < k + 1; i++)
    {
        h.x[i] = (c * b.x[i]) % N;
    }
    // g = v2o(h);
    printpol(h);
    printf(" =h in oinv2\n");
    return h;
}

vec vadd(vec a, vec b)
{
    int i;
    vec c = {0};

    // printf("deg=%d %d\n",deg(a),deg(b));

    for (i = 0; i < DEG; i++)
        c.x[i] = (a.x[i] + b.x[i]) % N;

    return c;
}

vec lsft(vec a)
{
    vec b = {0};
    int o = deg(a);

    for (int i = 0; i < o + 1; i++)
    {
        b.x[i + 1] = a.x[i];
    }
    // b.x[K*2]=0;

    return b;
}

vec rsft(vec a)
{
    vec b = {0};
    int o = deg(a);

    for (int i = 0; i < o + 1; i++)
        b.x[i] = a.x[i + 1];
    // b.x[0]=0;

    return b;
}

int mul = 0, mul2 = 0;
vec vmul(vec a, vec b)
{
    int i, j, k, l;
    vec c = {0};

    k = deg(a);
    l = deg(b);

    i = 0;
    while (i < k + 1)
    {
        for (j = 0; j < l + 1; j++)
        {
            if (a.x[i] > 0)
                c.x[i + j] = (c.x[i + j] + a.x[i] * b.x[j]) % N;
        }
        i++;
    }

    return c;
}

unsigned short vb[K * 2][N] = {0};
unsigned short gt[K * 2][K * 2] = {0};

void van(int kk)
{
    int i, j;

    printf("van der\n");

    for (i = 0; i < N; i++)
    {
        mat[i][0] = vb[0][i] = 1;
        printf("%d,", vb[0][i]);
    }
    printf("\n");

    // #pragma omp parallel for private(i, j)
    for (i = 1; i < kk; i++)
    {
        for (j = 0; j < N; j++)
        {
            vb[i][j] = mltn(i, j);
            printf("g%d,", vb[i][j]);
            mat[j][i] = vb[i][j];
        }
        printf("\n");
    }
}

void ogt(unsigned short pp[], int kk)
{
    int i, j;

    // #pragma omp parallel for private(i, j)
    for (i = 0; i < kk; i++)
    {
        for (j = 0; j < kk - i; j++)
        {
            gt[i][j + i] = g[j];
        }
    }
    for (i = 0; i < kk; i++)
    {
        for (j = 0; j < kk; j++)
            printf("h%d,", gt[i][j]);
        printf("\n");
    }
    // exit(1);
}

// 配列の値を係数として多項式に設定する
OP setpol(unsigned short f[], int n)
{
    OP g;
    vec v = {0};
    int i;

    for (i = 0; i < n; i++)
    {
        v.x[n - 1 - i] = f[i];
    }

    g = v2o(v);

    return g;
}

vec mkpol2(int s){
    vec a={0};
    int i;
    for(i=0;i<s;i++)
    a.x[i]=rand()%N;
a.x[s]=1;

return a;
}

OP mkpol()
{
    int i, j, k, flg, ii = 0;
    OP w = {0};

    do
    {
        // fail = 0;
        j = 0;
        k = 0;
        flg = 0;
        // l = 0;
        memset(g, 0, sizeof(g));
        // memset(ta, 0, sizeof(ta));
        memset(w.t, 0, sizeof(w));
        ginit();
        ii++;
        if (ii > 100)
        {
            printf("erro=%d\n", ii);
            exit(1);
        }

        for (i = 0; i < K; i++)
        {
            if (g[K - 1] > 0)
                flg = 1;
            if (i % 2 == 1 && g[i] > 0 && i < K)
                k++;
        }

        // 偶数項だけにならないようにする
        if ((k > 0 && flg == 0) || (k > 1 && flg == 1))
        // if(k>0)
        {
            w = setpol(g, K + 1);
            j = 1;
            // if(isquad(w)==-1)
            // exit(1);
        }
        // exit(1);

    } while (j == 0);

    printpol(o2v(w));
    printf(" ==g\n");
    // exit(1);

    return w;
}

unsigned short
v2a(oterm a)
{
    int j;

    if (a.a == 0)
        return 0;

    // printf("aa=%d\n",a.a);
    for (j = 0; j < M; j++)
    {
        if (j == a.a && a.a > 0)
        {
            // printf("j==%d\n",j);
            return j - 1;
        }
    }
    return 0;
}

void printsage(vec a)
{
    int i, j;
    oterm b;

    printf("poly=");
    for (i = 0; i < DEG; i++)
    {
        if (a.x[i] > 0)
        {
            b.a = a.x[i];
            b.n = i;
            j = v2a(b);
            printf("%d*X**%d+", b.a, i); // for GF(2^m)
        }
    }
}

// 多項式の代入値
unsigned short
trace(OP f, unsigned short x)
{
    unsigned short u = 0;
    vec v = o2v(f);
    int d = deg((v)) + 1;

    for (int i = 0; i < d; i++)
    {
        if (v.x[i] > 0)
            u = (u + (v.x[i] * mltn(i, x))) % N;
    }

    return u;
}

// リーディグタームを抽出(default)
oterm vLT(vec f)
{
    int i;
    oterm t = {0};

    // k = deg (o2v (f));
    for (i = 0; i < DEG; i++)
    {
        // printf("a=%d %d\n",f.t[i].a,f.t[i].n);
        if (f.x[i] > 0)
        {
            t.n = i;
            t.a = f.x[i];
        }
    }

    return t;
}

short inv(short a, short n)
{
    short d = n;
    short x = 0;
    short s = 1;
    while (a != 0)
    {
        short q = d / a;
        short r = d % a;
        d = a;
        a = r;
        short t = x - q * s;
        x = s;
        s = t;
    }
    short gcd = d; // $\gcd(a, n)$

    return ((x + n) % (n / d));
}

// aに何をかけたらbになるか
unsigned short
equ(unsigned short a, unsigned short b)
{
    // for(short i=0;i<N;i++)
    if (b == 0)
        return 0;
    if (a == 1)
        return b;

    return (inv(a, N) * b) % N;
}

// 多項式を単行式で割る
oterm vLTdiv(vec f, oterm t)
{
    oterm tt = {0}, s = {
                        0};

    tt = vLT(f);
    if (tt.n < t.n)
    {
        s.n = 0;
        s.a = 0;
    }
    else if (tt.n == t.n)
    {
        s.n = 0;
        s.a = equ(t.a, tt.a);
    }
    else if (tt.n > t.n)
    {
        s.n = tt.n - t.n;
        s.a = equ(t.a, tt.a);
        // printf("%u\n",s.a);
    }
    else if (t.n == 0 && t.a > 0)
    {
        s.a = (tt.a * inv(t.a, N)) % N;
        s.n = tt.n;
    }

    return s;
}

// 多項式を項ずつ掛ける
vec vterml(vec f, oterm t)
{
    // f = conv(f);
    // ssert(op_verify(f));
    int i;
    vec h = {0};

    // f=conv(f);
    // k = deg (o2v(f));

    for (i = 0; i < DEG; i++)
    {
        // h.t[i].n = f.t[i].n + t.n;
        if (f.x[i] > 0)
            h.x[i + t.n] = (f.x[i] * t.a) % N;
    }

    // h = conv(h);
    //  assert(op_verify(h));
    return h;
}

// 20200816:正規化したいところだがうまく行かない
// 多項式の足し算
vec vsub(vec a, vec b)
{
    vec c = {0};
    // int i, j, k, l = 0;
    vec h = {0}, f2 = {0}, g2 = {0};

    for (int i = 0; i < DEG; i++)
    {
        if (a.x[i] >= b.x[i])
            c.x[i] = (a.x[i] - b.x[i]) % N;
        if (a.x[i] < b.x[i])
            c.x[i] = (N + a.x[i] - b.x[i]) % N;
    }

    return c;
}

int vm = 0;
// 多項式の剰余を取る
vec vmod(vec f, vec g)
{
    vec h = {0};
    oterm b = {0}, c = {0};

    if (deg(g) == 0)
        return g;
    vm++;
    // printf("vmod-bl=%d k=%d\n",deg(f),deg(g));
    if (vLT(f).n < vLT(g).n)
    {
        //    exit(1);
        return f;
    }

    b = vLT(g);

    // printpol(f);
    // printf(" ==f\n");
    while (1)
    {
        // printf("@\n");
        c = vLTdiv(f, b);
        h = vterml(g, c);
        f = vsub(f, h);
        // printsage(g);
        if (deg((f)) == 0 || deg((h)) == 0)
        {
            break;
        }

        if (c.n == 0)
            break;
    }
    // printf("vmod-baka== %d %d\n",deg(f),deg(g));
    return f;
}

// int mul = 0, mul2 = 0;
vec vmul_2(vec a, vec b)
{
    int i, j, k, l;
    vec c = {0};
    if (deg(a) > 128 && deg(b) > 128)
        mul++;
    mul2++;

    k = deg(a);
    l = deg(b);

    for (i = 0; i < k + 1; i++)
    {
        for (j = 0; j < l + 1; j++)
        // if (a.x[i] > 0)
        {
            c.x[i + j] += (a.x[i] * b.x[j]) % N;
            c.x[i + j] %= N;
        }
    }

    return c;
}

// 多項式のべき乗
vec opow(vec f, int n)
{
    // int i;
    vec g = {0};

    g = f;

    for (int i = 1; i < n; i++)
        g = vmul(g, f);

    return g;
}

vec vpowmod(vec f, vec mod, int n)
{
    vec v = {0};
    vec ret = {0};

    v.x[0] = 1;
    ret = (v);
    while (n > 0)
    {
        if (n % 2 == 1)
            ret = vmod(vmul(ret, f), mod); // n の最下位bitが 1 ならば x^(2^i) をかける
        f = vmod(vmul(f, f), mod);
        n >>= 1; // n を1bit 左にずらす
    }
    return ret;
}

// gcd
vec ogcd(vec xx, vec yy)
{
    vec tt = {0}, tmp, h = {0};
    // ee.x[K] = 1;

    h.x[0] = 1;
    // h.x[0] = 0;
    if (deg((xx)) < deg((yy)))
    {
        tmp = xx;
        xx = yy;
        yy = tmp;
    }
    // tt = vmod(xx, yy);
    tt = vmod(xx, yy);
    while (deg(tt) > 0)
    {
        // printf("Oh!\n");
        xx = yy;
        yy = tt;
        if (deg(yy) > 0)
        {
            tt = vmod(xx, yy);
        }
        if (vLT(tt).a == 0)
            return yy;
    }
    if (vLT(yy).a == 0)
    {
        return tt;
    }
    else
    {
        return h;
    }
    //  return yy;
}

short diag(MTX a, int n)
{
    return (a.x[n][n] * a.x[n + 1][n + 1] - a.x[n][n + 1] * a.x[n + 1][n]) % N;
}

int resl(vec f, vec g)
{
    MTX a = {0};
    short dia[N] = {0};
    /*
    f.x[0]=16;
    f.x[1]=0;
    f.x[2]=4;
    f.x[3]=4;
    f.x[4]=1;
    g.x[0]=8;
    g.x[1]=9;
    g.x[2]=10;
    g.x[3]=9;
    printf("\n");
    */
    int n = deg(f), m = deg(g);
    if (n < m)
    {
        for (int i = 0; i < n + 1; i++)
        {
            for (int j = 0; j < m + 1; j++)
            {
                a.x[i + j][i] = f.x[n - j];
            }
        }
        for (int i = 0; i < n + m; i++)
        {
            for (int j = 0; j < n + m; j++)
            {
                a.x[j + i][i + m] = g.x[m - j];
            }
        }
    }
    if (n >= m)
    {
        for (int i = 0; i < m + 1; i++)
        {
            for (int j = 0; j < n + 1; j++)
            {
                a.x[i + j][i] = f.x[n - j];
            }
        }
        for (int i = 0; i < n + m + 1; i++)
        {
            for (int j = 0; j < n + m + 1; j++)
            {
                a.x[j + i][i + m] = g.x[m - j];
            }
        }
    }
    /*
    for(int i=0;i<n+m;i++){
        for(int j=0;j<m+n;j++)
        printf("%d,",a.x[i][j]);
        printf("\n");
    }
    printf("\n");
    */
    short tmp[N] = {0};
    int i, j, k, t;
    for (i = 0; i < m + n - 1; i++)
    {

        for (k = i; k < m + n - 1; k++)
        { // m+n
            // printf("%d ",k);
            t = a.x[k + 1][i];
            for (int j = i; j < n + m; j++)
            {
                tmp[j] = a.x[k + 1][j] - (a.x[i][j] * equ(a.x[i][i], a.x[k + 1][i])) % N;
                // printf("i=%d (j=%d k+1=%d) n=%d ks=%d %d %d t=%d =%d\n",i,j,k+1,a.x[k+1][j],(a.x[i][j]*equ(a.x[i][i],a.x[k+1][i]))%N,a.x[k][j],(a.x[i][j]),t,(N+tmp[j])%N);
            }
            // printf("\n");
            for (int j = 0; j < n + m; j++)
            {
                a.x[k + 1][j] = tmp[j];
                if (a.x[k + 1][j] < 0)
                    a.x[k + 1][j] = N + a.x[k + 1][j];
            }
            /*
            for(int u=0;u<n+m;u++){
                for(int v=0;v<n+m;v++)
                printf("%d ",a.x[u][v]);
                printf("\n");
            }
            printf(" %d %d %d\n",k,m+n,i);
            */
        }
        dia[i] = a.x[i][i];
    }

    int y = diag(a, n + m - 2);

    for (i = 0; i < m + n - 2; i++)
    {
        y = (y * dia[i]) % N;
        if (dia[i] == 0)
            return 0;
    }
    printf("y=%d\n", y);
    // exit(1);
    /*
    vec c=ogcd(f,g);
    if((deg(c)>0 && y>0)){ //} || (deg(c)==0 && y==0)){
    printsage(c);
    printf(" ==baka\n");
    printsage(f);
    printf(" ==f\n");
    printsage(g);
    printf(" ==g\n");
    exit(1);
    }
    */
    if (y > 0)
        return 0;
    if (y == 0)
        return -1;

    return 0;
}

int cnty = 0;
vec vpp(vec f, vec mod, int n)
{
    int i;
    vec s = {0};
    // t = f;
    s = f;
    printf("@\n");
    // 繰り返し２乗法
    for (i = 1; i < n; i++)
    {
        s = vmod(vmul_2(s, f), mod);
    }

    return s;
}

// GCD for decode
vec vgcd(vec xx, vec yy)
{
    vec tt;

    while (deg(yy) > 0)
    {
        tt = vmod(xx, yy);
        xx = yy;
        yy = tt;
    }
    if (yy.x[0] > 0)
        tt = kof2(yy.x[0], xx);
    printpol((yy));
    printf(" =========yy\n");
    printpol((tt));
    printf(" =========tt\n");

    return tt;
}

// 行列の逆行列を計算する関数
MTX inverseMatrix(MTX A, MTX A_inv, int start_row, int end_row)
{
    int i, j, k;
    short temp;

    // 単位行列を初期化
    for (i = 0; i < K / 2; i++)
    {
        for (j = 0; j < K / 2 + 1; j++)
        {
            A_inv.x[i][j] = (i == j) ? 1 : 0;
        }
    }

    // ガウス・ジョルダン法による逆行列の計算
    for (k = start_row; k < end_row; k++)
    {
        temp = A.x[k][k];
        for (j = 0; j < K / 2 + 1; j++)
        {
            A.x[k][j] = A.x[k][j] * oinv(temp, N) % N;
            A_inv.x[k][j] = A_inv.x[k][j] * oinv(temp, N) % N;
        }
        for (i = start_row; i < end_row; i++)
        {
            if (i != k)
            {
                temp = A.x[i][k] % N;
                for (j = 0; j < K / 2 + 1; j++)
                {
                    A.x[i][j] -= (A.x[k][j] * temp) % N;
                    A_inv.x[i][j] -= (A_inv.x[k][j] * temp) % N;
                }
            }
        }
    }
    vec x = {0};
    for (i = 0; i < K / 2; i++)
    {
        if (N > A.x[i][K / 2])
        {
            x.x[K / 2 - i] = (N - A.x[i][K / 2]) % N;
        }
        else
        {
            x.x[K / 2 - i] = A.x[i][K / 2] % N;
        }
    }
    for (int i = 0; i < K / 2; i++)
    {
        printf("in inverse ");
        for (int j = 0; j < K / 2; j++)
        {
            if (A_inv.x[i][j] < 0)
                A_inv.x[i][j] = N + A_inv.x[i][j]%N;
            printf("%d ", A_inv.x[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    // exit(1);

    /*
        x.x[0] = 1;

        vec vv = {0};
        OP pol = {0};
        pol = setpol(x.x, K / 2 + 1);
        printpol(o2v(pol));
        printf(" ==key\n");
        int key=0;
        for (i = 0; i < N; i++)
        {
            // v.x[i] = 0;
            if (trace(pol, i) % N == 0)
            {
                printf("error position=%d\n", i);
                vv.x[key++] = i;
            }
        }
        for (i = 0; i < K / 2; i++)
        {
            for (j = 0; j < K / 2 + 1; j++)
                printf("%d,", A_inv.x[i][j]%N);
            printf("\n");
        }
        // exit(1);
      */
    return A_inv;
}

// #define NN 16
vec sol(MTX a, int start, int end)
{
    int p, d;
    int i, j, k;
    vec v = {0};

    for (i = start; i < end; i++)
    {
        p = a.x[i][i];

        for (j = 0; j < (K / 2 + 1); j++)
        {
            a.x[i][j] = (a.x[i][j] * inv(p, N)) % N;
        }

        for (j = 0; j < K / 2; j++)
        {
            if (i != j)
            {
                d = a.x[j][i];

                for (k = i; k < (K / 2 + 1); k++)
                {
                    if (a.x[j][k] > (d * a.x[i][k]) % N)
                    {
                        a.x[j][k] -= (d * a.x[i][k]) % N;
                    }
                    else
                    {
                        a.x[j][k] = (N + (a.x[j][k] - (d * a.x[i][k]) % N)) % N;
                    }
                }
            }
        }
    }
    vec x = {0};
    for (i = start; i < end; i++)
    {
        if (N > a.x[i][K / 2])
        {
            x.x[K / 2 - i] = (N - a.x[i][K / 2]) % N;
        }
        else
        {
            x.x[K / 2 - i] = a.x[i][K / 2] % N;
        }
    }

    x.x[0] = 1;

    vec vv = {0};
    OP pol = {0};
    pol = setpol(x.x, K / 2 + 1);
    printpol(o2v(pol));
    printf(" ==key\n");
    int key = 0;
    for (i = 0; i < N; i++)
    {
        // v.x[i] = 0;
        if (trace(pol, i) % N == 0)
        {
            printf("error position=%d\n", i);
            vv.x[key++] = i;
        }
    }

    return vv;
}

// 多項式のべき乗余
vec opowmod(vec f, vec mod, int n)
{
    // int i, j = 0;
    vec g = f;
    printsage(mod);
    printf(" ma\n");
    // 繰り返し２乗法
    for (int i = 1; i < n; i++)
    {
        // f = vmul(f, f);
        g = vmul(g, f);
        if (deg(g) > deg(mod))
        {
            // printsage(g);
            // printf(" tadaima!\n");
            g = vmod(g, mod);
            // printsage(g);
            // printf(" tadaima2!\n");
        }
    }
    printsage(g);
    printf(" ==ge!\n");
    // exit(1);
    return g;
}

int is_equ(vec a, vec b)
{
    for (int i = 0; i < N * N; i++)
        if (a.x[i] != b.x[i])
            return -1;

    return 0;
}

// GF(2^m) then set m in this function.
int ben_or(vec f)
{
    int n; //, pid;

    vec s = {0}, u = {0}, r = {0};
    vec v = {0}; //, ff=o2v(f);
    // if GF(8192) is 2^m and m==13 or if GF(4096) and m==12 if GF(16384) is testing
    // int m = E;
    //  m=12 as a for GF(4096)=2^12 defined @ gloal.h or here,for example m=4 and GF(16)

    v.x[1] = 1;
    s = (v);
    // for (int i = 0; i < K / 2; i++)
    r = s;
    n = deg((f));

    if (vLT(f).n == 0)
    {
        printf("f==0\n");
        exit(1);
    }
    if (n == 0)
        return -1;

    // r(x)^{q^i} square pow mod
    for (int i = 0; i < K / 2; i++)
    {
        printf(":i=%d", i);
        // irreducible over GH(8192) 2^13
        // if(r.x[0]==65535)
        // return -1;
        // printsage(r);
        // printf(" --p\n");

        memset(r.x, 0, sizeof(r.x));
        v = vpowmod(v, f, N);
        r = v;
        // r.x[l]=1;

        u = vsub(r, (s));
        u = vmod(u, f);

        if (deg(u) > 0)
        {
            // printsage(u);
            // printf(" you\n");
            // printsage(f);
            printf(" me\n");
            u = ogcd(f, u);
            // int le=resl(f,u);
            // if(le==0 && deg(u)==0){
            //     printf("baka^^\n");
            // exit(1);
            // return -1;
            // }
            printf("you\n");
        }
        else
        {
            return -1;
        }
        if (deg(u) > 0) //  || vLT(u).a > 0)
        {
            // if(fequ(u,f)==1)
            {
                // flg[i]= -1;
                printf("ae\n");
                return -1;
            }
        }
    }

    return 0;
}

vec mkd(OP w, int kk, int start, int end)
{
    int i, j, k, l, ii = 0;

    unsigned short tr[N] = {0};
    unsigned short ta[N] = {0};
    vec v = {0}, pp = {0}, tt = {0};
    unsigned short po[K + 1] = {1, 0, 1, 0, 5};
    // vec w={0};
    vec r = {0};

aa:

    // printf("\n");
    memset(mat, 0, sizeof(mat));
    // 既約性判定のためのBen-Orアルゴリズム。拡大体にも対応している。デフォルトでGF(8192)
    // 既約多項式しか使わない。

    l = 0;
    ii = 0;
    // irreducible gvecpa code (既役多項式が必要なら、ここのコメントを外すこと。)

    w = mkpol();
    
    l = ben_or(o2v(w));
    while (l == -1)
        goto aa;
    printsage(o2v(w));
    printf("\n");
    //exit(1);
    
    //     printf("wwwwwww\n");
    //  exit(1);
    //  separable gvecpa code
    //  w = mkpol();
    r = o2v(w);
    //  r=vmul(w,w);
    memset(ta, 0, sizeof(ta));
    // w = setpol(g, K + 1);
    printpol((r));
    printf(" =poly\n");
    // exit(1);

    // 多項式の値が0でないことを確認
    for (int i = start; i < end; i++)
    {
        ta[i] = trace(w, i);
        if (ta[i] == 0)
        {
            printf("eval 0 @ %d\n", i);
            // fail = 1;
            // exit(1);
            goto aa;
        }
    }
    for (int i = start; i < end; i++)
    {
        tr[i] = inv(ta[i], N);
        // printf("%d,", tr[i]);
    }
    memset(g, 0, sizeof(g));
    // g[0] = 1;

    // 多項式を固定したい場合コメントアウトする。
    printpol(r);
    printf("\n");
    printsage((r));
    printf("\n");
    printf("sagemath で既約性を検査してください！\n");
    memset(v.x, 0, sizeof(v.x));
    //  v=rev(w);
    van(kk);
    //  v=(w);
    ogt(r.x, kk);
    // exit(1);
    //  wait();

    // #pragma omp parallel for

    printf("\nすげ、オレもうイキそ・・・\n");
    // keygen(g);
    // exit(1);

    for (int j = start; j < end; j++)
    {
        for (int i = 0; i < M; i++)
        {
            ma[i][j] = (vb[j][i] * tr[i]) % N;
        }
    }

    for (int i = start; i < end; i++)
    {
        for (int j = 0; j < M; j++)
        {
            for (int k = 0; k < K; k++)
            {
                mat[j][i] = (mat[j][i] + (gt[k][i] * ma[j][k])) % N;
            }
            //printf("c%d,", mat[j][i]);
        }
        //printf("\n");
    }

    /*
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < kk; j++)
            {
                mat[j][i] = vb[j][i];
            }
        }
    */
    // printf("\n");
    // exit(1);
    /*
    for( int j = 0; j < N; j++)
    {
        for( int i= 0; i < kk; i++)
            printf("%d,", mat[j][i]);
        printf("\n");
    }
    //exit(1);
    //wait();
*/

    return o2v(w);
}

void vv(int kk)
{
    int i, j;
    OP r = mkpol();
    unsigned short tr[N];
    unsigned short ta[N] = {0};

    printf("van der\n");

    for (i = 0; i < kk; i++)
    {
        for (j = 0; j < N; j++)
        {
            vb[i][j] = mltn(i, j);
        }
        // printf("\n");
    }

    int l = -1;
    vec pp = {0}, tt = {0};

aa:
    // exit(1);
    r = mkpol();

    for (i = 0; i < N; i++)
    {
        ta[i] = trace(r, i);
        if (ta[i] == 0)
        {
            printf("trace 0 @ %d\n", i);
            // fail = 1;
            goto aa;
        }
    }

    for (i = 0; i < N; i++)
    {
        tr[i] = inv(ta[i], N);
        // printf("%d,", tr[i]);
    }

    printf("\nすげ、オレもうイキそ・・・\n");
    // keygen(g);
    // exit(1);

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < kk; j++)
        {
            mat[i][j] = (vb[j][i] * tr[i]) % N;
        }
    }
    for (i = 0; i < K; i++)
    {
        for (j = 0; j < N; j++)
            printf("c%d,", mat[j][i]);
        printf("\n");
    }
}

void mkerr(unsigned short *z1, int num)
{
    int j, l;

    j = 0;

    memset(z1, 0, sizeof(2 * N));

    while (j < num)
    {
        l = rand() % (N - 1);
        // printf ("l=%d\n", l);
        if (0 == z1[l] && l > 0)
        {
            z1[l] = 2; //rand()%N;
            // printf("l=%d\n", l);
            if(z1[l]>0)
            j++;
        }
    }
}

OP synd(unsigned short zz[], int kk)
{
    unsigned short syn[K] = {0}, s = 0;
    int i, j;
    OP f = {0};

    printf("in synd2\n");

    for (i = 0; i < kk; i++)
    {
        syn[i] = 0;
        s = 0;
        // #pragma omp parallel num_threads(16)
        for (j = 0; j < N; j++)
        {
            s = (s + (zz[j] * mat[j][i])) % N;
        }
        syn[i] = s;
        // printf ("syn%d,", syn[i]);
    }
    // printf ("\n");

    f = setpol(syn, kk);
    printpol(o2v(f));
    printf(" syn============= %d\n", deg(o2v(f)));
    //  exit(1);

    return f;
}

// chen探索
vec chen(vec f)
{
    vec e = {0};
    int i, n, x = 0, count = 0;
    unsigned short z;

    n = deg((f));
    for (x = 0; x < N; x++)
    {
        z = 0;
        for (i = 0; i < n + 1; i++)
        {
            if (f.x[i] > 0)
                z += (mltn(i, x) * f.x[i]) % N;
        }
        if (z % N == 0)
        {
            e.x[count] = x;
            count++;
            printf("change %d\n", (x));
        }
    }

    return e;
}

typedef struct
{
    vec f;
    vec g;
    vec h;
} ymo;

vec pmul(vec a, vec b)
{
    int i, j, k, l;
    vec c = {0};

    k = deg(a) + 1;
    l = deg(b) + 1;
    printf("k=%d,l=%d", k, l);
    // exit(1);
    for (int i = 0; i < k; i++)
    {
        for (int j = 0; j < l; j++)
            if (a.x[i] > 0)
            {
                c.x[i + j] = (c.x[i + j] + a.x[i] * b.x[j]) % N;
                // printf("%d=c ",c.x[i+j]);
            }
        // printf("\n");
    }
    /*
    printf("\n");
    printpol(v2o(c));
    printf(" ==c\n");
    printpol(v2o(a));
    printf(" ==a\n");
    printpol(v2o(b));
    printf(" ==b\n");
    // exit(1);
    */
    return c;
}

vec bms(unsigned short s[])
{
    int L = 0, m = -1, d[K] = {0}, k = 0, i, e;
    vec f = {0}, g = {0}, h, v;

    f.x[0] = g.x[0] = 1;

    while (k <= (2 * T - 1))
    {
        e = 0;
        for (i = 0; i < L; i++)
            e = (e+f.x[i]*s[k - i])%N;

        d[k] = (f.x[i]*s[k - i] + e)%N; // s[k] ^ e;
        if (d[k] > 0)
        {
            h = f;
            memset(v.x, 0, sizeof(v.x));
            v.x[k - m] = 1;

            unsigned short a;
            a = (m < 0) ? 1 : oinv(d[m],N);
            f = vadd(f, vmul(kof2((d[k]* a), g), v));
            if (L <= k / 2)
            {
                L = k + 1 - L;
                m = k;
                g = h;
            }
        }
        k++;
    }

    return f;
}

ymo bm_itr(unsigned short s[])
{
    vec U1[2][2] = {0}, U2[2][2][2] = {0}, null = {0};
    int i, j, k;
    ymo t = {0};

    U2[0][0][0].x[0] = 1;       // f[0];
    U2[0][0][1].x[0] = 0;       // fai[0];
    U2[0][1][0].x[0] = 0;       // g[0];
    U2[0][1][1].x[0] = N - (1); // thi[0];
    int m = 0, d = 0, p = (2 * d - m - 1)%N, myu = 0;
    printf("m=%d d=%d myu=%d p=%d\n", m, d, myu, p);
    for (m = 0; m < K; m++)
    {
        d = deg(U2[0][0][0]);
        p = (2 * d - m - 1)%N;
        myu = 0;
        for (int i = 0; i <= d; i++)
            myu = (myu + U2[0][0][0].x[i] * s[i + (m - d)]) % N;

        printf("m=%d ad=%d myu=%d p=%d\n", m, d, myu, p);
        memset(U1, 0, sizeof(U1));
        if (myu == 0 || p >= 0)
        {
            U1[0][0].x[0] = 1;
            U1[0][1].x[p] = N - (myu);
            U1[1][0].x[0] = 0;
            U1[1][1].x[0] = 1;
            // exit(1);
        }
        else if (myu > 0 && p < 0)
        {
            if (p < 0)
            {
                p = -1 * (p);
            }
            U1[0][0].x[p] = 1;
            U1[0][1].x[0] = N - (myu);
            U1[1][0].x[0] = oinv(myu, N);
            U1[1][1].x[0] = 0;
        }
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                for (int k = 0; k < 2; k++)
                    U2[1][i][j] = (vadd((U2[1][i][j]), (pmul(U1[i][k], U2[0][k][j]))));
            }
        }
        memcpy(U2[0], U2[1], sizeof(U2[0]));
        memset(U2[1], 0, sizeof(U2[1]));
    }
    t.f = U2[0][0][0];
    t.g = U2[0][1][0];
    t.h = U2[0][0][1];
    if (deg(t.f) == T)
    {
        printsage((t.f));
        printf(" ==chen00\n");
        return t;
    }
    else
    {
        t.f = U2[1][0][0];
        printsage((t.f));
        printf("baka\n");
        exit(1);
    }
}

// 行列の掛け算関数
void matrix_multiply(short A[MATRIX_SIZE][MATRIX_SIZE], short B[MATRIX_SIZE][MATRIX_SIZE], short *C, int start_row, int end_row)
{
    for (int i = start_row; i < end_row; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            short sum = 0.0;
            for (int k = 0; k < MATRIX_SIZE; k++)
            {
                sum += A[i][k] * B[k][j];
            }
            C[i * MATRIX_SIZE + j] = sum;
        }
    }
}

int ipow(int b, int n)
{
    int l = 1;

    if (n == 0)
        return 1;

    for (int i = 0; i < n; i++)
        l = b * l % N;

    return l;
}

MTX ver(MTX x)
{
    MTX d = {0};
    int det = 0;

    det = oinv(((x.x[0][0] * x.x[1][1] - x.x[0][1] * x.x[1][0])), N);
    printf("det=%d %d %d\n", det, (x.x[0][0] * x.x[1][1] - x.x[0][1] * x.x[1][0]), oinv(-6, N));
    d.x[0][0] = (det * x.x[1][1]) % N;
    d.x[1][1] = x.x[0][0] * det % N;
    d.x[0][1] = (N - x.x[0][1] * det % N);
    d.x[1][0] = (N - x.x[1][0] * det % N);
    printf("X=%d %d\n", x.x[0][0], x.x[0][1]);
    printf("X=%d %d\n", x.x[1][0], x.x[1][1]);
    printf("A=%d %d\n", d.x[0][0], d.x[0][1]);
    printf("A=%d %d\n", d.x[1][0], d.x[1][1]);
    // exit(1);

    return d;
}

vec recv(MTX t, vec v)
{
    vec x = {0};
    int i, j, k;
    for (i = 0; i < K/2; i++)
    {
        for (k = 0; k < K/2; k++)
            x.x[i] += (v.x[k] * t.x[i][k])%N;
        x.x[i] %= N;
        if(x.x[i]<0)
        x.x[i]=(x.x[i]+N)%N;
    }

    return x;
}

MTX mmul(MTX s, MTX t)
{
    MTX y = {0};

    y.x[0][0] = (s.x[0][0] * t.x[0][0] + s.x[0][1] * t.x[1][0]) % N;
    y.x[0][1] = (s.x[0][0] * t.x[0][1] + s.x[0][1] * t.x[1][1]) % N;
    y.x[1][0] = (s.x[1][0] * t.x[0][0] + s.x[1][1] * t.x[1][0]) % N;
    y.x[1][1] = (s.x[1][0] * t.x[0][1] + s.x[1][1] * t.x[1][1]) % N;
    printf("%d %d\n", y.x[0][0], y.x[0][1]);
    printf("%d %d\n", y.x[1][0], y.x[1][1]);

    return y;
}

vec ev(vec x,vec v)
{
    int i, j, k;
    MTX mmk = {0};
    MTX inv_A = {0};
    vec tx = {0};

    for (i = 0; i < K / 2; i++)
    {
        for (int j = 0; j < K / 2; j++)
        {
            mmk.x[j][i] = mat[x.x[i]][j];
            printf("%d %df", mat[x.x[j]][i], x.x[j]);
        }
        printf("\n");
    }
    printf("\n(");
    for (i = 0; i < K / 2; i++)
        printf("%d ", x.x[i]);
    printf(")\n");
    // exit(1);
    //mmk.x[0][K / 2] = 2;
    //mmk.x[1][K / 2] = 5;
    for (i = 0; i < K / 2; i++)
    {
        mmk.x[i][K / 2] = v.x[i];
        for (int j = 0; j < K / 2; j++)
            printf("%d^ ", mmk.x[i][j]);
        printf("\n");
    }

    // tx.x[1]=v.x[1];
    for (i = 0; i < K / 2; i++)
        tx.x[i] = v.x[i];
        //tx.x[1] = 5; //v.x[i];
    // inv_A=ver(mmk);
    inv_A = inverseMatrix(mmk, inv_A, 0, K / 2);
    printf("inv %d %d %d\n", inv_A.x[0][0], inv_A.x[0][1], inv_A.x[0][2]);
    printf("inv %d %d %d\n", inv_A.x[1][0], inv_A.x[1][1], inv_A.x[1][2]);
    tx = recv(inv_A, tx);
    for (i = 0; i < K / 2; i++)
        printf("error value is %d\n", tx.x[i]);

        return tx;
}

//モニック多項式にする
vec coeff(vec f, unsigned short d)
{
  int i, j, k;
  vec a, b;

  k = deg((f)) + 1;
  for (i = 0; i < k; i++)
    f.x[i] = (f.x[i]*oinv(d,N))%N;

  return f;
}


//多項式の商を取る
vec vdiv(vec f, vec g)
{

  int i = 0, j, n, k;
  vec h = {0}, e = {0}, tt = {0};
  oterm a, b = {0}, c = {0};

  if (vLT(f).n == 0 && vLT(g).a == 0)
  {
    printf("baka^\n");
    //return f;
    exit(1);
  }
  if (vLT(g).a == 0)
  {
    exit(1);
  }
  if (vLT(g).n == 0 && vLT(g).a > 1)
    return coeff(f, vLT(g).a);

  k = deg(g);
  b = vLT(g);
  if (b.a == 1 && b.n == 0)
    return f;
  if (b.a == 0 && b.n == 0)
  {
    printf("baka in vdiv\n");
    exit(1);
  }
  if (deg((f)) < deg((g)))
  {
    return f;
    //  a=LT(f);
  }

  i = 0;
  while (vLT(f).n > 0 && vLT(g).n > 0)
  {
    c = vLTdiv(f, b);
    assert(c.n < DEG);
    tt.x[c.n] = c.a;
    //i++;

    h = vterml(g, c);

    f = vsub(f, h);
    if (deg((f)) == 0 || deg((g)) == 0)
    {
      //printf ("blake2\n");
      break;
    }

    if (c.n == 0)
      break;
  }

  // tt は逆順に入ってるので入れ替える
  return tt;
}


// invert of polynomial
vec vinv(vec a, vec n)
{
  vec d = {0}, x = {0}, s = {0}, q = {0}, r = {0}, t = {0}, u = {0}, v = {0}, w = {0}, tt = {0}, gcd = {0}, tmp = {0};
  oterm b = {0};
  vec vv = {0}, xx = {0},aa={0},bb={0},cc={0};

  if (deg((a)) > deg((n)))
  {
    tmp = a;
    a = n;
    n = tmp;
    printf("baka_i\n");
    //exit (1);
  }
  if (vLT(a).a == 0)
  {
    printf(" a ga 0\n");
    exit(1);
  }

  tt = n;

  d = n;
  x.x[0] = 0;
  s.x[0] = 1;
  while (deg((a)) > 0)
  {
    if (deg((a)) > 0)
      r = vmod(d, a);
    if (vLT(a).a == 0)
      break;
    if (vLT(a).a > 0)
      q = vdiv(d, a);

    d = a;
    a = r;
    t = vsub(x, vmul_2(q, s));

    x = s;
    s = t;
  }
  d = a;
  a = r;

  x = s;
  s = t;

  gcd = d; // $\gcd(a, n)$
  printpol((gcd));
  printf(" =========gcd\n");

  printf("before vinv\n");
  //w=tt;

  b = vLT(w);

  aa=(x);
  bb=(n);
  cc = vadd(aa, bb);
  v=(cc);
  w = tt;
  if (vLT(v).n > 0 && vLT(w).n > 0)
  {
    u = vmod(v, w);
  }
  else
  {
    //printpol (o2v (v));
    printf(" v===========\n");
    //printpol (o2v (x));
    printf(" x==0?\n");
    //printpol (o2v (n));
    printf(" n==0?\n");

    exit(1);
  }
  //caution !!
  if (vLT(u).a > 0 && vLT(d).a > 0)
  {
    u = vdiv(u, d);
  }

  if (vLT(u).a == 0 || vLT(d).a == 0)
  {
    printf("inv div u or d==0\n");
    // exit(1);
  }
  if (vLT(u).a == 0)
  {
    printf("no return at u==0\n");
    exit(1);
  }

  return u;
}

int main()
{
    int i, u = 0;
    unsigned short s[K + 1] = {0}, z1[N] = {0};
    
    srand(clock());

    vec v = {0}, x = mkpol2(15),r1=mkpol2(12),r2=mkpol2(12),m={1,2,3,4,5,6,7,8},I={1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1};
    OP f = {0},y={0};
    printf("%d %d %d\n", 3, oinv(3, N), 3 * oinv(3, N) % N);
    int le = 1;
    for (i = 1; i < 65; i++)
    {
        for (int j = 0; j < 13; j++)
            le = (13 * le) % 128;
        printf("le=%d %d\n", le, i);
    }
    //v=vmod(vadd(vmul_2(x,r1),m),x);
    v=vmod(vmul(vinv(m,I),m),I);
    //v=vmod(c, (x));
    printpol(v);
    printf("\n");
    //exit(1);
    // mkg(K); // Goppa Code (EEA type)
    // van(K); // RS-Code generate
    // mkd(f, K);
    // vv(K);           // Goppa Code's Parity Check (Berlekamp type)

    // resl(v,x);
    // exit(1);

    mkd(f, K, 0, K);

    while (1)
    {
        // for(i=0;i<T;i++)
        // z1[i]=2;
        memset(z1, 0, sizeof(z1));
         mkerr(z1, T);    // generate error vector
        // for (int i = 0; i < T; i++)
        //    z1[i] = i+1;
        //z1[1] = 1;
        //z1[2] = 2;
        //z1[3] = 3;
        //z1[4] = 1;
        f = synd(z1, K); // calc syndrome
        x = o2v(f);      // transorm to vec
        vec r={0};
        //vec r = bms(x.x);    // Berlekamp-Massey Algorithm
        for(i=0;i<K;i++)
        v.x[K-i-1]=x.x[i];
        ymo y=bm_itr(v.x);
        //for(i=0;i<T;i++)
        //r.x[K-1-i]=y.f.x[i];
        x=chen(y.f);
        //chen(r);
        // exit(1);
         for(i=0;i<N;i++)
         if(z1[i]>0)
         printf("i=%d\n",i);
        //exit(1);
         // mkd(1);
        MTX b = {0};

        //for (i = 0; i < K; i++)
        //    v.x[K - 1 - i] = x.x[i];
        //printpol((v));
        //printf(" ==synpol\n");
        printpol((v));
        printf(" ==synpol\n");
        //exit(1);
/*
        for (i = 0; i < K / 2; i++)
        {
            for (int j = 0; j < K / 2 + 1; j++)
            {
                b.x[i][j] = v.x[i + j];
                // printf("%d,",b.x[i][i+j]);
            }
            // printf("\n");
        }
        printf("\n");
        for (i = 0; i < K / 2; i++)
        {
            for (int j = 0; j < K / 2 + 1; j++)
                printf("e%d,", b.x[i][j]);
            printf("\n");
        }
*/
        //x = sol(b, 0, K / 2);
        /*
        for (i = 0; i < N; i++)
        {
            if (z1[i] > 0 && x.x[i] == 0)
            {
                printf("baka=%d %d %d\n", i, z1[i], x.x[i]);
                // exit(1);
            }
        }
        */
        x=ev(x,v);
        //exit(1);
        int flg = 0,yo=0;
        for (i = 0; i < N; i++)
        {
            if (z1[i] > 0)
            {
                printf("(correcting ,original) = (%d, %d) %d\n", x.x[yo], z1[i], i);
                yo++;
                flg++;
            }
        }

        if (flg == T)
            break;

        break;
    }
    return 0;
}
