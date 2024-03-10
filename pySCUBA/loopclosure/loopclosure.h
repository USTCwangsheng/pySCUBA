/*
 * loopclosure.h
 *
 *  Created on: 2016年12月7日
 *      Author: hyliu
 */

#ifndef LOOPCLOSURE_LOOPCLOSURE_H_
#define LOOPCLOSURE_LOOPCLOSURE_H_


#define PRINT_LEVEL 0
#define	MAX_ORDER 16
#define MAXPOW 32
#define	SMALL_ENOUGH 1.0e-18

namespace NSPloopclosure {

typedef  	struct	p {
  int	ord;
  double	coef[MAX_ORDER+1];
} poly;

void initialize_sturm(double *tol_secant, int *max_iter_sturm, int *max_iter_secant);
void solve_sturm(const int *p_order, int *n_root, double *poly_coeffs, double *roots);
double hyper_tan(double a, double x);
static int modp(poly *u, poly *v, poly *r);
int buildsturm(int ord, poly *sseq);
int numroots(int np, poly *sseq, int *atneg, int *atpos);
int numchanges(int np, poly *sseq, double a);
void sbisect(int np, poly *sseq, double min, double max, int atmin, int atmax, double *roots);
double evalpoly(int ord, double *coef, double x);
int modrf(int ord, double *coef, double	a, double b, double *val);


double dot_product(double va[], double vb[]);
void matmul(double ma[3][3], double mb[3], double mc[3]);
double sign(double a, double b);
void solve_3pep_poly(double r_n1[], double r_a1[], double r_a3[], double r_c3[], double r_soln_n[][3][3], double r_soln_a[][3][3], double r_soln_c[][3][3], int *n_soln);
void initialize_loop_closure(double b_len[], double b_ang[], double t_ang[]);
void get_input_angles(int *n_soln, double r_n1[], double r_a1[], double r_a3[], double r_c3[]);
void test_two_cone_existence_soln(double tt, double kx, double et, double ap, int *n_soln, char cone_type[]);
void get_poly_coeff(double poly_coeff[]);
void poly_mul_sub2(double u1[][5], double u2[][5], double u3[][5], double u4[][5], int p1[], int p2[], int p3[], int p4[], double u5[][5], int p5[]);
void poly_mul2(double u1[][5], double u2[][5], int p1[], int p2[], double u3[][5], int p3[]);
void poly_sub2(double u1[][5], double u2[][5], int p1[], int p2[], double u3[][5], int p3[]);
void poly_mul_sub1(double u1[], double u2[], double u3[], double u4[], int p1, int p2, int p3, int p4, double u5[], int *p5);
void poly_mul1(double u1[], double u2[], int p1, int p2, double u3[], int *p3);
void poly_sub1(double u1[], double u2[], int p1, int p2, double u3[], int *p3);
void coord_from_poly_roots(int *n_soln, double roots[], double r_n1[], double r_a1[], double r_a3[], double r_c3[], double r_soln_n[][3][3], double r_soln_a[][3][3], double r_soln_c[3][3][3]);
double calc_t2(double t0);
double calc_t1(double t0, double t2);
void calc_dih_ang(double r1[], double r2[], double r3[], double *angle);
void calc_bnd_ang(double r1[], double r2[], double *angle);
void cross(double p[], double q[], double s[]);
void quaternion(double axis[], double quarter_ang, double p[]);
void rotation_matrix(double q[], double U[3][3]);
}




#endif /* LOOPCLOSURE_LOOPCLOSURE_H_ */
