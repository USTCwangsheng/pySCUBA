/*
 * tripepclosure.cpp
 *	a wrap-up of the analytical loop closure code of
 *	 Chaok Seok, Evangelos Coutsias, Matthew Jacobson, and Ken Dill
 *  Created on: 2016年12月7日
 *      Author: hyliu
 */
#include "loopclosure/loopclosure.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
namespace NSPloopclosure {
/*
double dot_product(double va[], double vb[]);
void matmul(double ma[][], double mb[], double mc[]);
double sign(double a, double b);
void solve_3pep_poly(double r_n1[], double r_a1[], double r_a3[], double r_c3[], double r_soln_n[][][], double r_soln_a[][][], double r_soln_c[][][], int *n_soln);
void initialize_loop_closure(double b_len[], double b_ang[], double t_ang[]);
void get_input_angles(int *n_soln, double r_n1[], double r_a1[], double r_a3[], double r_c3[]);
void test_two_cone_existence_soln(double tt, double kx, double et, double ap, int *n_soln, char cone_type[]);
void get_poly_coeff(double poly_coeff[]);
void poly_mul_sub2(double u1[][], double u2[][], double u3[][], double u4[][], int p1[], int p2[], int p3[], int p4[], double u5[][], int p5[]);
void poly_mul2(double u1[][], double u2[][], int p1[], int p2[], double u3[][], int p3[]);
void poly_sub2(double u1[][], double u2[][], int p1[], int p2[], double u3[][], int p3[]);
void poly_mul_sub1(double u1[], double u2[], double u3[], double u4[], int p1, int p2, int p3, int p4, double u5[], int *p5);
void poly_mul1(double u1[], double u2[], int p1, int p2, double u3[], int *p3);
void poly_sub1(double u1[], double u2[], int p1, int p2, double u3[], int *p3);
void coord_from_poly_roots(int *n_soln, double roots[], double r_n1[], double r_a1[], double r_a3[], double r_c3[], double r_soln_n[][][], double r_soln_a[][][], double r_soln_c[][][]);
double calc_t2(double t0);
double calc_t1(double t0, double t2);
void calc_dih_ang(double r1[], double r2[], double r3[], double *angle);
void calc_bnd_ang(double r1[], double r2[], double *angle);
void cross(double p[], double q[], double s[]);
void quaternion(double axis[], double quarter_ang, double p[]);
void rotation_matrix(double q[], double U[][]);
*/
//!----------------------------------------------------------------------
//! Copyright (C) 2003 
//!      Chaok Seok, Evangelos Coutsias, Matthew Jacobson, and Ken Dill
//!      UCSF and Univeristy of New Mexico
//! Witten by Chaok Seok 2003.  
//!----------------------------------------------------------------------------
//!----------------------------------------------------------------------------
//!*************** Tripeptide Loop Closure Algorithm *****************
//! files to be compiled with:
//!    tripeptide_closure.f90
//!    sturm.c
//! 
//!*******************************************************************
//! subroutine  initialize_loop_closure(b_len, b_ang, t_ang)
//!*******************************************************************
//! connectivity of atoms:
//!   N1-A1-C1-N2-A2-C2-N3-A3-C3
//!
//! input:
//!
//  * b_len(1:6): bond lengths (A1-C1, C1-N2, ..., N3-A3)
//!  * b_ang(1:7): bond angles (N1-A1-C1, A1-C1-N2, ..., N3-A3-C3)
//!  * t_ang(1:2): torsion angles (A1-C1-N2-A2, A2-C2-N3-A3)
//!*******************************************************************
//!
//!*******************************************************************
//! subroutine solv_3pep_poly(r_n1, r_a1, r_a3, r_c3, &
//!     r_soln_n, r_soln_a, r_soln_c, n_soln)
//!*******************************************************************
//! input: 
//!  * r_n1(3), r_a1(3), r_a3(3), r_c3(3): 
//!       Cartesian coordinates of N and CA atoms of the first residue and
//!        CA and C atoms of the last (third) residue.
//! output:
//!  * n_soln: number of alternative loop closure solutions.
//!  * r_soln_n(3,3,8), r_soln_a(3,3,8), r_soln_c(3,3,8): 
//!       Cartesian coordinates of loop closure solutions. 
//!       first dimension: x, y, z component
//!       second dim: residue number
//!       third dim: solution number
//!*******************************************************************
//!----------------------------------------------------------------------------
//MODULE tripep_closure
//!----------------------------------------------------------------------------
//  integer, parameter :: dp = kind(1.0d0)
//  real(dp), parameter :: pi = 3.141592653589793238462643383279502884197d0
  #define  pi 3.141592653589793238462643383279502884197e0
//  real(dp), parameter :: two_pi=2.0d0*pi, deg2rad = pi/180.0d0, rad2deg = 180.0d0/pi
  #define two_pi 2.0e0*pi
  #define deg2rad pi/180.0e0
  #define rad2deg 180.0e0/pi
  #define max(a,b) ((a) > (b))? (a) : (b)
  #define min(a,b) ((a) < (b))? (a) : (b)
//  integer, parameter :: max_soln = 16
  const int max_soln = 16;
//  integer, parameter :: deg_pol = 16
  const int deg_pol = 16;
//  integer, parameter :: print_level = 0
//  int print_level = 1;
  int print_level = 0;
//  ! parameters for tripeptide loop (including bond lengths & angles)
//  real(dp) :: len0(6), b_ang0(7), t_ang0(2)
  double len0[6], b_ang0[7], t_ang0[2];
//  real(dp) :: aa13_min_sqr, aa13_max_sqr
  double aa13_min_sqr, aa13_max_sqr;
//  real(dp) :: delta(0:3), xi(3), eta(3), alpha(3), theta(3)
  double delta[4], xi[3], eta[3], alpha[3], theta[3];
//  real(dp) :: cos_alpha(3), sin_alpha(3), cos_theta(3), sin_theta(3)
  double cos_alpha[3], sin_alpha[3], cos_theta[3], sin_theta[3];
//  real(dp) :: cos_delta(0:3), sin_delta(0:3)
  double cos_delta[4], sin_delta[4];
//  real(dp) :: cos_xi(3), cos_eta(3), sin_xi(3), sin_eta(3)
  double cos_xi[3], cos_eta[3], sin_xi[3], sin_eta[3];
//  real(dp) :: r_a1a3(3), r_a1n1(3), r_a3c3(3)
  double r_a1a3[3], r_a1n1[3], r_a3c3[3];
//  real(dp) :: b_a1a3(3), b_a1n1(3), b_a3c3(3)
  double b_a1a3[3], b_a1n1[3], b_a3c3[3];
//  real(dp) :: len_na(3), len_ac(3), len_aa(3)
  double len_na[3], len_ac[3], len_aa[3];
//  ! used for polynomial coefficients
//  real(dp) :: C0(0:2,3), C1(0:2,3), C2(0:2,3)
  double C0[3][3], C1[3][3], C2[3][3];
//  real(dp) :: Q(0:16,0:4), R(0:16,0:2)
  double Q[5][17], R[3][17];
//CONTAINS

double dot_product(double va[3], double vb[3])
{
 return va[0]*vb[0] + va[1]*vb[1] + va[2]*vb[2];
}
  
void matmul(double ma[3][3], double mb[3], double mc[3])
 {
  int i, j;

  for(i=0;i<3;i++)
   {
    mc[i] = 0.;
    for(j=0;j<3;j++)
      mc[i] += ma[i][j]*mb[j];
   }
  return;
 }

double sign(double a, double b)
 {
  if (b>=0.)
    return fabs(a);
  else
    return -fabs(a);
 }
//!-----------------------------------------------------------------------
//subroutine solv_3pep_poly(r_n1, r_a1, r_a3, r_c3, &
//     r_soln_n, r_soln_a, r_soln_c, n_soln)
void solve_3pep_poly(double r_n1[3], double r_a1[3], double r_a3[3], double r_c3[3], double r_soln_n[max_soln][3][3], double r_soln_a[max_soln][3][3], double r_soln_c[max_soln][3][3], int *n_soln)
 {
//  implicit none
//  real(dp), intent(in) :: r_n1(3), r_a1(3), r_a3(3), r_c3(3)
//  integer, intent(out) :: n_soln
//  real(dp), intent(out) :: r_soln_n(:,:,:), r_soln_a(:,:,:), r_soln_c(:,:,:)
//  real(dp) :: poly_coeff(0:deg_pol), roots(max_soln)
  double poly_coeff[deg_pol+1], roots[max_soln];

//  call get_input_angles(n_soln, r_n1, r_a1, r_a3, r_c3)
  get_input_angles(n_soln, r_n1, r_a1, r_a3, r_c3);

//  if (n_soln == 0) then
//     return
//  end if
  if (*n_soln == 0)
    return;

//  call get_poly_coeff(poly_coeff)
  get_poly_coeff(poly_coeff);

//  call solve_sturm(deg_pol, n_soln, poly_coeff, roots)
  solve_sturm(&deg_pol, n_soln, poly_coeff, roots);

//  if (n_soln == 0) then
//!     print*, 'return 2'
//     return
//  end if
  if (*n_soln == 0) return;
//    exit(1);

//  call coord_from_poly_roots(n_soln, roots, r_n1, r_a1, r_a3, r_c3, r_soln_n, r_soln_a, r_soln_c)
  coord_from_poly_roots(n_soln, roots, r_n1, r_a1, r_a3, r_c3, r_soln_n, r_soln_a, r_soln_c);

  return;
//end subroutine solv_3pep_poly
 }
//!-----------------------------------------------------------------------
//subroutine initialize_loop_closure(b_len, b_ang, t_ang)
void initialize_loop_closure(double b_len[6], double b_ang[7], double t_ang[2])
 {
//!-----------------------------------------------------------------------
//! Input angles for the given bond lengths and angles 
//!-----------------------------------------------------------------------
//  implicit none
//  real(dp), intent(in) :: b_len(6), b_ang(7), t_ang(2)
//  real(dp) :: len1, len2, a_min, a_max
  double len1, len2, a_min, a_max;
//  real(dp), dimension(3) :: axis, rr_a1, rr_c1, rr_n2, rr_a2, rr_n2a2_ref, rr_c1a1
  double axis[3], rr_a1[3], rr_c1[3], rr_n2[3], rr_a2[3], rr_n2a2_ref[3], rr_c1a1[3];
//  real(dp), dimension(3) :: rr_a1a2, dr, bb_c1a1, bb_a1a2, bb_a2n2
  double rr_a1a2[3], dr[3], bb_c1a1[3], bb_a1a2[3], bb_a2n2[3];
//  real(dp), dimension(3) :: r_a1, r_n1
//  double r_a1[3], r_n1[3];
//  real(dp) :: p(4), Us(3,3)
  double p[4], Us[3][3];
  double mulpro[3];
  double tmp_val[3];
//  real(dp), parameter :: tol_secant = 1.0d-15
  double tol_secant = 1.0e-15;
//  integer, parameter :: max_iter_sturm = 100, max_iter_secant = 20
  int max_iter_sturm = 100;
  int max_iter_secant = 20;
//  integer :: i
  int i, j;
  
//  call initialize_sturm(tol_secant, max_iter_sturm, max_iter_secant)
  initialize_sturm(&tol_secant, &max_iter_sturm, &max_iter_secant);

//  len0(1:6) = b_len(1:6)
  for(i=0;i<6;i++)
    len0[i] = b_len[i];
//  b_ang0(1:7) = b_ang(1:7)
  for(i=0;i<7;i++)
    b_ang0[i] = b_ang[i];
//  t_ang0(1:2) = t_ang(1:2)
  for(i=0;i<2;i++)
    t_ang0[i] = t_ang[i];

//  rr_c1(1:3) = 0.0d0
  for(i=0;i<3;i++)
    rr_c1[i] = 0.;
//  axis(1:3) = (/ 1.0d0, 0.0d0, 0.0d0 /)
  axis[0] = 1.;
  axis[1] = 0.;
  axis[2] = 0.;

//  do i = 0, 1
  for(i=0;i<2;i++)
   {
//     rr_a1(1:3) = (/ cos(b_ang0(3*i+2))*len0(3*i+1), sin(b_ang0(3*i+2))*len0(3*i+1), 0.0d0 /)
     rr_a1[0] = cos(b_ang0[3*i+1])*len0[3*i];
     rr_a1[1] = sin(b_ang0[3*i+1])*len0[3*i];
     rr_a1[2] = 0.0e0;
//     rr_n2(1:3) = (/ len0(3*i+2), 0.0d0, 0.0d0 /)
     rr_n2[0] = len0[3*i+1];
     rr_n2[1] = 0.0e0;
     rr_n2[2] = 0.0e0;
//     rr_c1a1(:) = rr_a1(:) - rr_c1(:)
     for(j=0;j<3;j++)
       rr_c1a1[j] = rr_a1[j] - rr_c1[j];
//     rr_n2a2_ref(1:3) = (/ -cos(b_ang0(3*i+3))*len0(3*i+3), sin(b_ang0(3*i+3))*len0(3*i+3), 0.0d0 /)
     rr_n2a2_ref[0] = -cos(b_ang0[3*i+2])*len0[3*i+2];
     rr_n2a2_ref[1] = sin(b_ang0[3*i+2])*len0[3*i+2];
     rr_n2a2_ref[2] = 0.0e0;
//     call quaternion(axis, t_ang0(i+1)*0.25d0, p)
     quaternion(axis, t_ang0[i]*0.25e0, p);
//     call rotation_matrix(p, Us)
     rotation_matrix(p, Us);
//     rr_a2(:) =  matmul(Us, rr_n2a2_ref) + rr_n2(:)
//     rr_a1a2(:) = rr_a2(:) - rr_a1(:)
//     dr(:) = rr_a1a2(:)
     matmul(Us, rr_n2a2_ref, mulpro);
     for(j=0;j<3;j++)
      {
       rr_a2[j] =  mulpro[j] + rr_n2[j];
       rr_a1a2[j] = rr_a2[j] - rr_a1[j];
       dr[j] = rr_a1a2[j];
      }
//     len2 = dot_product(dr, dr)
     len2 = dot_product(dr, dr);
//     len1 = sqrt(len2)
     len1 = sqrt(len2);
//     ! len_aa
//     len_aa(i+2) = len1
     len_aa[i+1] = len1;
//     bb_c1a1(:) = rr_c1a1(:)/len0(3*i+1)
//     bb_a1a2(:) = rr_a1a2(:)/len1
//     bb_a2n2(:) = (rr_n2(:) - rr_a2(:))/len0(3*i+3)
     for(j=0;j<3;j++)
      {
       bb_c1a1[j] = rr_c1a1[j]/len0[3*i];
       bb_a1a2[j] = rr_a1a2[j]/len1;
       bb_a2n2[j] = (rr_n2[j] - rr_a2[j])/len0[3*i+2];
      }
//     ! xi
//     call calc_bnd_ang(-bb_a1a2, bb_a2n2, xi(i+2))
     for(j=0;j<3;j++)
       tmp_val[j] = -bb_a1a2[j];
     calc_bnd_ang(tmp_val, bb_a2n2, &xi[i+1]);
//     ! eta
//     call calc_bnd_ang(bb_a1a2, -bb_c1a1, eta(i+1))
     for(j=0;j<3;j++)
       tmp_val[j] = -bb_c1a1[j];
     calc_bnd_ang(bb_a1a2, tmp_val, &eta[i]);
//     ! delta: pi -  dih of N(1)CA(1)CA(3)C(3)
//     call calc_dih_ang(bb_c1a1, bb_a1a2, bb_a2n2, delta(i+1))
     calc_dih_ang(bb_c1a1, bb_a1a2, bb_a2n2, &delta[i+1]);
//     delta(i+1) = pi - delta(i+1)
     delta[i+1] = pi - delta[i+1];
   }
//  end do

//  a_min = b_ang(4) - (xi(2) + eta(2))
  a_min = b_ang[3] - (xi[1] + eta[1]);
//  a_max = min(b_ang(4) + (xi(2) + eta(2)), pi)
  a_max = min((b_ang[3] + (xi[1] + eta[1])), pi);

//  ! min/max of base length
//!  print*, 'len1, len3=', len_aa(2:3)
//  printf("len1, len3= %9.5f %9.5f\n", len_aa[1], len_aa[2]);
//!  print*, 'a_min, a_max=', a_min*rad2deg, a_max*rad2deg
//  printf("a_min, a_max= %9.5f %9.5f\n", a_min*rad2deg, a_max*rad2deg);
//  aa13_min_sqr = len_aa(2)**2 + len_aa(3)**2 - 2.0d0*len_aa(2)*len_aa(3)*cos(a_min)
  aa13_min_sqr = pow(len_aa[1],2) + pow(len_aa[2],2) - 2.0e0*len_aa[1]*len_aa[2]*cos(a_min);
//  aa13_max_sqr = len_aa(2)**2 + len_aa(3)**2 - 2.0d0*len_aa(2)*len_aa(3)*cos(a_max)
  aa13_max_sqr = pow(len_aa[1],2) + pow(len_aa[2],2) - 2.0e0*len_aa[1]*len_aa[2]*cos(a_max);
//!  print*, 'aa13_min_sqr,aa13_max_sqr', aa13_min_sqr,aa13_max_sqr
//  printf("aa13_min_sqr,aa13_max_sqr %9.5f %9.5f\n", aa13_min_sqr, aa13_max_sqr);

//end subroutine initialize_loop_closure
 }
//!-----------------------------------------------------------------------
//subroutine get_input_angles(n_soln, r_n1, r_a1, r_a3, r_c3)
void get_input_angles(int *n_soln, double r_n1[3], double r_a1[3], double r_a3[3], double r_c3[3])
 {
//!-----------------------------------------------------------------------
//! Input angles and vectors (later used in coordinates) 
//!-----------------------------------------------------------------------
//  implicit none
//  real(dp), intent(in) :: r_n1(:), r_a1(:), r_a3(:), r_c3(:)
//  integer, intent(out) :: n_soln
//  real(dp) :: dr_sqr
  double dr_sqr;
  double tmp_val[3];
//  integer :: i
  int i;
//  character(len=2) :: cone_type
  char cone_type[2];
//!-----------------------------------------------------------------------

//  n_soln = max_soln
  *n_soln = max_soln;

//  ! vertual bond
//  r_a1a3(:) = r_a3(:) - r_a1(:) 
  for(i=0;i<3;i++)
    r_a1a3[i] = r_a3[i] - r_a1[i]; 
//  dr_sqr = dot_product(r_a1a3,r_a1a3)
  dr_sqr = dot_product(r_a1a3,r_a1a3);
//  len_aa(1) = sqrt(dr_sqr)
  len_aa[0] = sqrt(dr_sqr);

//  if (dr_sqr < aa13_min_sqr .or. dr_sqr > aa13_max_sqr) then
//     n_soln = 0
//!     print*, 'return 0'
//!     print*, sqrt(dr_sqr), sqrt(aa13_min_sqr), sqrt(aa13_max_sqr)
//     return
//  end if
  if ((dr_sqr < aa13_min_sqr) || (dr_sqr > aa13_max_sqr))
   {
    *n_soln = 0;
    return;
   }

//  ! bond lengths
//  r_a1n1(:) = r_n1(:) - r_a1(:)
  for(i=0;i<3;i++)
    r_a1n1[i] = r_n1[i] - r_a1[i];
//  len_na(1) = sqrt(dot_product(r_a1n1,r_a1n1))
  len_na[0] = sqrt(dot_product(r_a1n1,r_a1n1));
//  len_na(2) = len0(3)
  len_na[1] = len0[2];
//  len_na(3) = len0(6)
  len_na[2] = len0[5];
//  r_a3c3(:) = r_c3(:) - r_a3(:)
  for(i=0;i<3;i++)
    r_a3c3[i] = r_c3[i] - r_a3[i];
//  len_ac(1) = len0(1)
  len_ac[0] = len0[0];
//  len_ac(2) = len0(4)
  len_ac[1] = len0[3];
//  len_ac(3) = sqrt(dot_product(r_a3c3,r_a3c3))
  len_ac[2] = sqrt(dot_product(r_a3c3,r_a3c3));

//  ! unit vectors
//  b_a1n1(:) = r_a1n1(:)/len_na(1)
//  b_a3c3(:) = r_a3c3(:)/len_ac(3)
//  b_a1a3(:) = r_a1a3(:)/len_aa(1)
  for(i=0;i<3;i++)
   {
    b_a1n1[i] = r_a1n1[i]/len_na[0];
    b_a3c3[i] = r_a3c3[i]/len_ac[2];
    b_a1a3[i] = r_a1a3[i]/len_aa[0];
   }

//  ! delta(3): dih of N(1)CA(1)CA(3)C(3)
//  call calc_dih_ang(-b_a1n1, b_a1a3, b_a3c3, delta(3))
  for(i=0;i<3;i++)
    tmp_val[i] = -b_a1n1[i];
  calc_dih_ang(tmp_val, b_a1a3, b_a3c3, &delta[3]);
//  delta(0) = delta(3)
  delta[0] = delta[3];

//  ! xi(1) 
//  call calc_bnd_ang(-b_a1a3, b_a1n1, xi(1))
  for(i=0;i<3;i++)
    tmp_val[i] = -b_a1a3[i];
  calc_bnd_ang(tmp_val, b_a1n1, &xi[0]);
 
//  ! eta(3)
//  call calc_bnd_ang(b_a1a3, b_a3c3, eta(3))
  calc_bnd_ang(b_a1a3, b_a3c3, &eta[2]);

//  do i = 1, 3
  for(i=0;i<3;i++)
   {
//     cos_delta(i) = cos(delta(i))
     cos_delta[i+1] = cos(delta[i+1]);
//     sin_delta(i) = sin(delta(i))
     sin_delta[i+1] = sin(delta[i+1]);
//     cos_xi(i) = cos(xi(i))
//     cos_xi(i) = cos(xi(i))
     cos_xi[i] = cos(xi[i]);
//     sin_xi(i) = sin(xi(i))
     sin_xi[i] = sin(xi[i]);
//     sin_xi(i) = sin(xi(i))
     sin_xi[i] = sin(xi[i]);
//     cos_eta(i) = cos(eta(i))
     cos_eta[i] = cos(eta[i]);
//     cos_eta(i) = cos(eta(i))
     cos_eta[i] = cos(eta[i]);
//     sin_eta(i) = sin(eta(i))
     sin_eta[i] = sin(eta[i]);
//     sin_eta(i) = sin(eta(i))
     sin_eta[i] = sin(eta[i]);
//  end do
   }
//  cos_delta(0) = cos_delta(3)
  cos_delta[0] = cos_delta[3];
//  sin_delta(0) = sin_delta(3)
  sin_delta[0] = sin_delta[3];

//  ! theta (N, CA, C) bond angle
//  theta(1) = b_ang0(1)
  theta[0] = b_ang0[0];
//  theta(2) = b_ang0(4)
  theta[1] = b_ang0[3];
//  theta(3) = b_ang0(7)
  theta[2] = b_ang0[6];
//  do i = 1, 3
//     cos_theta(i) = cos(theta(i))
//  end do
  for(i=0;i<3;i++)
    cos_theta[i] = cos(theta[i]);

//  ! alpha 
//  cos_alpha(1) = -(len_aa(1)**2 + len_aa(2)**2 - len_aa(3)**2)/(2.0d0*len_aa(1)*len_aa(2))
  cos_alpha[0] = -(pow(len_aa[0],2) + pow(len_aa[1],2) - pow(len_aa[2],2))/(2.0e0*len_aa[0]*len_aa[1]);
//  alpha(1) = acos(cos_alpha(1))
  alpha[0] = acos(cos_alpha[0]);
//  sin_alpha(1) = sin(alpha(1))
  sin_alpha[0] = sin(alpha[0]);
//  cos_alpha(2) = (len_aa(2)**2 + len_aa(3)**2 - len_aa(1)**2)/(2.0d0*len_aa(2)*len_aa(3))
  cos_alpha[1] = (pow(len_aa[1],2) + pow(len_aa[2],2) - pow(len_aa[0],2))/(2.0e0*len_aa[1]*len_aa[2]);
//  alpha(2) = acos(cos_alpha(2))
  alpha[1] = acos(cos_alpha[1]);
//  sin_alpha(2) = sin(alpha(2))
  sin_alpha[1] = sin(alpha[1]);
//  alpha(3) = pi - alpha(1) + alpha(2)
  alpha[2] = pi - alpha[0] + alpha[1];
//  cos_alpha(3) = cos(alpha(3))
  cos_alpha[2] = cos(alpha[2]);
//  sin_alpha(3) = sin(alpha(3))
  sin_alpha[2] = sin(alpha[2]);

//  if (print_level > 0) then
  if (print_level > 0)
   {
//     write(*,'(a,3f9.4)') 'xi = ', xi(1:3)*rad2deg
     printf("xi = %9.4f%9.4f%9.4f\n", xi[0]*rad2deg, xi[1]*rad2deg, xi[2]*rad2deg);
//     write(*,'(a,3f9.4)') 'eta = ', eta(1:3)*rad2deg
     printf("eta = %9.4f%9.4f%9.4f\n", eta[0]*rad2deg, eta[1]*rad2deg, eta[2]*rad2deg);
//     write(*,'(a,3f9.4)') 'delta = ', delta(1:3)*rad2deg
     printf("delta = %9.4f%9.4f%9.4f\n", delta[1]*rad2deg, delta[2]*rad2deg, delta[3]*rad2deg);
//     write(*,'(a,3f9.4)') 'theta = ', theta(1:3)*rad2deg
     printf("theta = %9.4f%9.4f%9.4f\n", theta[0]*rad2deg, theta[1]*rad2deg, theta[2]*rad2deg);
//     write(*,'(a,3f9.4)') 'alpha = ', alpha(1:3)*rad2deg
     printf("alpha = %9.4f%9.4f%9.4f\n", alpha[0]*rad2deg, alpha[1]*rad2deg, alpha[2]*rad2deg);
//  end if
   }

//  ! check for existence of soln
//  do i = 1, 3
  for(i=0;i<3;i++)
   {
//     call test_two_cone_existence_soln(theta(i), xi(i), eta(i), alpha(i), &
//          n_soln, cone_type)
    test_two_cone_existence_soln(theta[i], xi[i], eta[i], alpha[i], n_soln, cone_type);
//     if (n_soln == 0) then
//        print*, 'return 1', i
//        return
//     end if
    if (*n_soln == 0)
      return;
//  end do
   }

  return;
//end subroutine get_input_angles
 }
//!-----------------------------------------------------------------------
//subroutine test_two_cone_existence_soln(tt, kx, et, ap, n_soln, cone_type)
void test_two_cone_existence_soln(double tt, double kx, double et, double ap, int *n_soln, char cone_type[2])
 {
//  implicit none
//  real(dp), intent(in) :: tt, kx, et, ap
//  integer, intent(out) :: n_soln
//  character(len=2), intent(out) :: cone_type
//  character(len=2) :: case_type
//  real(dp) :: at, ex, abs_at, ap1, kx1, et1
  double at, ex, abs_at, ap1, kx1, et1;
//  real(dp) :: cos_tx1, cos_tx2, cos_te1, cos_te2, cos_ea1, cos_ea2, cos_xa1, cos_xa2
  double cos_tx1, cos_tx2, cos_te1, cos_te2, cos_ea1, cos_ea2, cos_xa1, cos_xa2;
//  logical :: s1, s2, t1, t2, complicated = .false.
  int s1, s2, t1, t2;
  int complicated = 0;
//  real(dp), parameter :: half_pi = 0.5d0*pi

//  n_soln = max_soln
  *n_soln = max_soln;
 
//  ap1 = ap
  ap1 = ap;
//  kx1 = kx
  kx1 = kx;
//  et1 = et
  et1 = et;
  
//  at = ap1 - tt
  at = ap1 - tt;
//  ex = kx1 + et1
  ex = kx1 + et1;
//  abs_at = abs(at)
  abs_at = fabs(at);

//  ! case of no soln
//  if (abs_at > ex) then
//     n_soln = 0
//     return
//  end if
  if (abs_at > ex)
   {
    *n_soln = 0;
    return;
   }

//  if (complicated) then
//     ! find type of intersection
//     cos_tx1 = cos(tt+kx1)
//     cos_tx2 = cos(tt-kx1)
//     cos_te1 = cos(tt+et1)
//     cos_te2 = cos(tt-et1)
//     cos_ea1 = cos(et1+ap1)
//     cos_ea2 = cos(et1-ap1)
//     cos_xa1 = cos(kx1+ap1)
//     cos_xa2 = cos(kx1-ap1)
//     s1 = .false.; s2 = .false.; t1 = .false.; t2 = .false. 
//     if ((cos_te1-cos_xa2)*(cos_te1-cos_xa1) <= 0.0d0) s1 = .true.
//     if ((cos_te2-cos_xa2)*(cos_te2-cos_xa1) <= 0.0d0) s2 = .true.
//     if ((cos_tx1-cos_ea2)*(cos_tx1-cos_ea1) <= 0.0d0) t1 = .true.
//     if ((cos_tx2-cos_ea2)*(cos_tx2-cos_ea1) <= 0.0d0) t2 = .true.
//  end if
  if (complicated)
   {
//    cos_tx1 = cos(tt+kx1)
    cos_tx1 = cos(tt+kx1);
//    cos_tx2 = cos(tt-kx1)
    cos_tx2 = cos(tt-kx1);
//    cos_te1 = cos(tt+et1)
    cos_te1 = cos(tt+et1);
//    cos_te2 = cos(tt-et1)
    cos_te2 = cos(tt-et1);
//    cos_ea1 = cos(et1+ap1)
    cos_ea1 = cos(et1+ap1);
//    cos_ea2 = cos(et1-ap1)
    cos_ea2 = cos(et1-ap1);
//    cos_xa1 = cos(kx1+ap1)
    cos_xa1 = cos(kx1+ap1);
//    cos_xa2 = cos(kx1-ap1)
    cos_xa2 = cos(kx1-ap1);
//    s1 = .false.; s2 = .false.; t1 = .false.; t2 = .false. 
    s1 = 0;
    s2 = 0;
    t1 = 0;
    t2 = 0;
//    if ((cos_te1-cos_xa2)*(cos_te1-cos_xa1) <= 0.0d0) s1 = .true.
    if ((cos_te1-cos_xa2)*(cos_te1-cos_xa1) <= 0.0e0)
      s1 = 0;
//    if ((cos_te2-cos_xa2)*(cos_te2-cos_xa1) <= 0.0d0) s2 = .true.
    if ((cos_te2-cos_xa2)*(cos_te2-cos_xa1) <= 0.0e0)
      s2 = 0;
//    if ((cos_tx1-cos_ea2)*(cos_tx1-cos_ea1) <= 0.0d0) t1 = .true.
    if ((cos_tx1-cos_ea2)*(cos_tx1-cos_ea1) <= 0.0e0)
      t1 = 0;
//    if ((cos_tx2-cos_ea2)*(cos_tx2-cos_ea1) <= 0.0d0) t2 = .true.
    if ((cos_tx2-cos_ea2)*(cos_tx2-cos_ea1) <= 0.0e0)
      t2 = 0;
   }

  return;
//end subroutine test_two_cone_existence_soln
 }
//!-----------------------------------------------------------------------
//subroutine get_poly_coeff(poly_coeff)
void get_poly_coeff(double poly_coeff[deg_pol+1])
 {
//  implicit none
//  real(dp), intent(out) :: poly_coeff(0:deg_pol)
//  integer :: i, j
  int i, j;
//  real(dp) :: A0, A1, A2, A3, A4, A21, A22, A31, A32, A41, A42
  double A0, A1, A2, A3, A4, A21, A22, A31, A32, A41, A42;
//  real(dp) :: B0(3), B1(3), B2(3), B3(3), B4(3), B5(3), B6(3), B7(3), B8(3)
  double B0[3], B1[3], B2[3], B3[3], B4[3], B5[3], B6[3], B7[3], B8[3];
//  real(dp), dimension(0:4,0:4) :: u11, u12, u13, u31, u32, u33
  double u11[5][5], u12[5][5], u13[5][5], u31[5][5], u32[5][5], u33[5][5];
//  real(dp), dimension(0:4,0:4) :: um1, um2, um3, um4, um5, um6, q_tmp
  double um1[5][5], um2[5][5], um3[5][5], um4[5][5], um5[5][5], um6[5][5], q_tmp[5][5];
//  integer, dimension(2) :: p1, p3, p_um1, p_um2, p_um3, p_um4, p_um5, p_um6, p_Q
  int p1 [2], p3[2], p_um1[2], p_um2[2], p_um3[2], p_um4[2], p_um5[2], p_um6[2], p_Q[2];
//  integer :: p2, p4, p_f1, p_f2, p_f3, p_f4, p_f5, p_f6, p_f7, &
//       p_f8, p_f9, p_f10, p_f11, p_f12, p_f13, p_f14, p_f15, p_f16, p_f17, &
//       p_f18, p_f19, p_f20, p_f21, p_f22, p_f23, p_f24, p_f25, p_f26, p_f27
  int p2, p4, p_f1, p_f2, p_f3, p_f4, p_f5, p_f6, p_f7, p_f8, p_f9;
  int p_f10, p_f11, p_f12, p_f13, p_f14, p_f15, p_f16, p_f17, p_f18;
  int p_f19, p_f20, p_f21, p_f22, p_f23, p_f24, p_f25, p_f26;
//  integer :: p_final
  int p_final;
//  real(dp), dimension(0:16) :: f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, &
//       f12, f13, f14, f15, f16, f17, f18, f19, f20, f21, f22, f23, f24, f25, f26, f27
  double f1[17], f2[17], f3[17], f4[17], f5[17], f6[17], f7[17], f8[17], f9[17];
  double f10[17], f11[17], f12[17], f13[17], f14[17], f15[17], f16[17], f17[17], f18[17];
  double f19[17], f20[17], f21[17], f22[17], f23[17], f24[17], f25[17], f26[17];

//  ! A0, B0
//  do i = 1, 3
  for(i=0;i<3;i++)
   {
//     A0 = cos_alpha(i)*cos_xi(i)*cos_eta(i) - cos_theta(i)
    A0 = cos_alpha[i]*cos_xi[i]*cos_eta[i] - cos_theta[i];
//     A1 = -sin_alpha(i)*cos_xi(i)*sin_eta(i)
    A1 = -sin_alpha[i]*cos_xi[i]*sin_eta[i];
//     A2 = sin_alpha(i)*sin_xi(i)*cos_eta(i)
    A2 = sin_alpha[i]*sin_xi[i]*cos_eta[i];
//     A3 = sin_xi(i)*sin_eta(i)
    A3 = sin_xi[i]*sin_eta[i];
//     A4 = A3*cos_alpha(i)
    A4 = A3*cos_alpha[i];
//     j = i - 1
    j = i;
//     A21 = A2*cos_delta(j)
    A21 = A2*cos_delta[j];
//     A22 = A2*sin_delta(j)
    A22 = A2*sin_delta[j];
//     A31 = A3*cos_delta(j)
    A31 = A3*cos_delta[j];
//     A32 = A3*sin_delta(j)
    A32 = A3*sin_delta[j];
//     A41 = A4*cos_delta(j)
    A41 = A4*cos_delta[j];
//     A42 = A4*sin_delta(j)
    A42 = A4*sin_delta[j];
//     B0(i) = A0 + A22 + A31
    B0[i] = A0 + A22 + A31;
//     B1(i) = 2.0d0*(A1 + A42)
    B1[i] = 2.0e0*(A1 + A42);
//     B2(i) = 2.0d0*(A32 - A21)
    B2[i] = 2.0e0*(A32 - A21);
//     B3(i) = -4.0d0*A41
    B3[i] = -4.0e0*A41;
//     B4(i) = A0 + A22 - A31
    B4[i] = A0 + A22 - A31;
//     B5(i) = A0 - A22 - A31
    B5[i] = A0 - A22 - A31;
//     B6(i) = -2.0d0*(A21 + A32)
    B6[i] = -2.0e0*(A21 + A32);
//     B7(i) = 2.0d0*(A1 - A42)
    B7[i] = 2.0e0*(A1 - A42);
//     B8(i) = A0 - A22 + A31
    B8[i] = A0 - A22 + A31;
//  end do
   }

//  ! C0i
//  i = 1
  i = 0;
//  C0(0:2,i) = (/ B0(i), B2(i), B5(i) /);
  C0[i][0] = B0[i];
  C0[i][1] = B2[i];
  C0[i][2] = B5[i];
//  C1(0:2,i) = (/ B1(i), B3(i), B7(i) /);
  C1[i][0] = B1[i];
  C1[i][1] = B3[i];
  C1[i][2] = B7[i];
//  C2(0:2,i) = (/ B4(i), B6(i), B8(i) /);
  C2[i][0] = B4[i];
  C2[i][1] = B6[i];
  C2[i][2] = B8[i];
//  do i = 2, 3
  for(i=1;i<3;i++)
   {
//     C0(0:2,i) = (/ B0(i), B1(i), B4(i) /)
    C0[i][0] = B0[i];
    C0[i][1] = B1[i];
    C0[i][2] = B4[i];
//     C1(0:2,i) = (/ B2(i), B3(i), B6(i) /)
    C1[i][0] = B2[i];
    C1[i][1] = B3[i];
    C1[i][2] = B6[i];
//     C2(0:2,i) = (/ B5(i), B7(i), B8(i) /)
    C2[i][0] = B5[i];
    C2[i][1] = B7[i];
    C2[i][2] = B8[i];
//  end do
   }

//  ! first determinant
//  do i = 0, 2
  for(i=0;i<3;i++)
   {
//     u11(i,0) = C0(i,1)
     u11[0][i] = C0[0][i];
//     u12(i,0) = C1(i,1)
     u12[0][i] = C1[0][i];
//     u13(i,0) = C2(i,1)
     u13[0][i] = C2[0][i];
//     u31(0,i) = C0(i,2) STRANGE !!!
     u31[i][0] = C0[1][i];
//     u32(0,i) = C1(i,2) STRANGE !!!
     u32[i][0] = C1[1][i];
//     u33(0,i) = C2(i,2) STRANGE !!!
     u33[i][0] = C2[1][i];
//  end do
   }

//  p1(1:2) = (/ 2, 0 /)
  p1[0] = 2;
  p1[1] = 0;
//  p3(1:2) = (/ 0, 2 /)
  p3[0] = 0;
  p3[1] = 2;

//  call poly_mul_sub2(u32, u32, u31, u33, p3, p3, p3, p3, um1, p_um1)
  poly_mul_sub2(u32, u32, u31, u33, p3, p3, p3, p3, um1, p_um1);
//  call poly_mul_sub2(u12, u32, u11, u33, p1, p3, p1, p3, um2, p_um2)
  poly_mul_sub2(u12, u32, u11, u33, p1, p3, p1, p3, um2, p_um2);
//  call poly_mul_sub2(u12, u33, u13, u32, p1, p3, p1, p3, um3, p_um3)
  poly_mul_sub2(u12, u33, u13, u32, p1, p3, p1, p3, um3, p_um3);
//  call poly_mul_sub2(u11, u33, u31, u13, p1, p3, p3, p1, um4, p_um4)
  poly_mul_sub2(u11, u33, u31, u13, p1, p3, p3, p1, um4, p_um4);
//  call poly_mul_sub2(u13, um1, u33, um2, p1, p_um1, p3, p_um2, um5, p_um5)
  poly_mul_sub2(u13, um1, u33, um2, p1, p_um1, p3, p_um2, um5, p_um5);
//  call poly_mul_sub2(u13, um4, u12, um3, p1, p_um4, p1, p_um3, um6, p_um6)
  poly_mul_sub2(u13, um4, u12, um3, p1, p_um4, p1, p_um3, um6, p_um6);
//  call poly_mul_sub2(u11, um5, u31, um6, p1, p_um5, p3, p_um6, q_tmp, p_Q)
  poly_mul_sub2(u11, um5, u31, um6, p1, p_um5, p3, p_um6, q_tmp, p_Q);

//  Q(0:4,0:4) = q_tmp(0:4,0:4)
  for(i=0;i<5;i++)
    for(j=0;j<5;j++)
      Q[i][j] = q_tmp[i][j];

//  ! second determinant
//  R(:,:) = 0.0d0
  for(i=0;i<3;i++)
    for(j=0;j<17;j++)
      R[i][j] = 0.;
//  R(0:2,0) = C0(0:2,3)
//  R(0:2,1) = C1(0:2,3)
//  R(0:2,2) = C2(0:2,3)
  for(i=0;i<3;i++)
   {
    R[0][i] = C0[2][i];
    R[1][i] = C1[2][i];
    R[2][i] = C2[2][i];
   }
//  p2 = 2
  p2 = 2;
//  p4 = 4
  p4 = 4;

//  call poly_mul_sub1(R(:,1), R(:,1), R(:,0), R(:,2), p2, p2, p2, p2, f1, p_f1)
  poly_mul_sub1(R[1], R[1], R[0], R[2], p2, p2, p2, p2, f1, &p_f1);
//  call poly_mul1(R(:,1), R(:,2), p2, p2, f2, p_f2)
  poly_mul1(R[1], R[2], p2, p2, f2, &p_f2);
//  call poly_mul_sub1(R(:,1), f1, R(:,0), f2, p2, p_f1, p2, p_f2, f3, p_f3)
  poly_mul_sub1(R[1], f1, R[0], f2, p2, p_f1, p2, p_f2, f3, &p_f3);
//  call poly_mul1(R(:,2), f1, p2, p_f1, f4, p_f4)
  poly_mul1(R[2], f1, p2, p_f1, f4, &p_f4);
//  call poly_mul_sub1(R(:,1), f3, R(:,0), f4, p2, p_f3, p2, p_f4, f5, p_f5)
  poly_mul_sub1(R[1], f3, R[0], f4, p2, p_f3, p2, p_f4, f5, &p_f5);

//  call poly_mul_sub1(Q(:,1), R(:,1), Q(:,0), R(:,2), p4, p2, p4, p2, f6, p_f6)
  poly_mul_sub1(Q[1], R[1], Q[0], R[2], p4, p2, p4, p2, f6, &p_f6);
//  call poly_mul_sub1(Q(:,2), f1, R(:,2), f6, p4, p_f1, p2, p_f6, f7, p_f7)
  poly_mul_sub1(Q[2], f1, R[2], f6, p4, p_f1, p2, p_f6, f7, &p_f7);
//  call poly_mul_sub1(Q(:,3), f3, R(:,2), f7, p4, p_f3, p2, p_f7, f8, p_f8)
  poly_mul_sub1(Q[3], f3, R[2], f7, p4, p_f3, p2, p_f7, f8, &p_f8);
//  call poly_mul_sub1(Q(:,4), f5, R(:,2), f8, p4, p_f5, p2, p_f8, f9, p_f9)
  poly_mul_sub1(Q[4], f5, R[2], f8, p4, p_f5, p2, p_f8, f9, &p_f9);

//  call poly_mul_sub1(Q(:,3), R(:,1), Q(:,4), R(:,0), p4, p2, p4, p2, f10, p_f10)
  poly_mul_sub1(Q[3], R[1], Q[4], R[0], p4, p2, p4, p2, f10, &p_f10);
//  call poly_mul_sub1(Q(:,2), f1, R(:,0), f10, p4, p_f1, p2, p_f10, f11, p_f11)
  poly_mul_sub1(Q[2], f1, R[0], f10, p4, p_f1, p2, p_f10, f11, &p_f11);
//  call poly_mul_sub1(Q(:,1), f3, R(:,0), f11, p4, p_f3, p2, p_f11, f12, p_f12)
  poly_mul_sub1(Q[1], f3, R[0], f11, p4, p_f3, p2, p_f11, f12, &p_f12);

//  call poly_mul_sub1(Q(:,2), R(:,1), Q(:,1), R(:,2), p4, p2, p4, p2, f13, p_f13)
  poly_mul_sub1(Q[2], R[1], Q[1], R[2], p4, p2, p4, p2, f13, &p_f13);
//  call poly_mul_sub1(Q(:,3), f1, R(:,2), f13, p4, p_f1, p2, p_f13, f14, p_f14)
  poly_mul_sub1(Q[3], f1, R[2], f13, p4, p_f1, p2, p_f13, f14, &p_f14);
//  call poly_mul_sub1(Q(:,3), R(:,1), Q(:,2), R(:,2), p4, p2, p4, p2, f15, p_f15)
  poly_mul_sub1(Q[3], R[1], Q[2], R[2], p4, p2, p4, p2, f15, &p_f15);
//  call poly_mul_sub1(Q(:,4), f1, R(:,2), f15, p4, p_f1, p2, p_f15, f16, p_f16)
  poly_mul_sub1(Q[4], f1, R[2], f15, p4, p_f1, p2, p_f15, f16, &p_f16);
//  call poly_mul_sub1(Q(:,1), f14, Q(:,0), f16, p4, p_f14, p4, p_f16, f17, p_f17)
  poly_mul_sub1(Q[1], f14, Q[0], f16, p4, p_f14, p4, p_f16, f17, &p_f17);

//  call poly_mul_sub1(Q(:,2), R(:,2), Q(:,3), R(:,1), p4, p2, p4, p2, f18, p_f18)
  poly_mul_sub1(Q[2], R[2], Q[3], R[1], p4, p2, p4, p2, f18, &p_f18);
//  call poly_mul_sub1(Q(:,1), R(:,2), Q(:,3), R(:,0), p4, p2, p4, p2, f19, p_f19)
  poly_mul_sub1(Q[1], R[2], Q[3], R[0], p4, p2, p4, p2, f19, &p_f19);
//  call poly_mul_sub1(Q(:,3), f19, Q(:,2), f18, p4, p_f19, p4, p_f18, f20, p_f20)
  poly_mul_sub1(Q[3], f19, Q[2], f18, p4, p_f19, p4, p_f18, f20, &p_f20);
//  call poly_mul_sub1(Q(:,1), R(:,1), Q(:,2), R(:,0), p4, p2, p4, p2, f21, p_f21)
  poly_mul_sub1(Q[1], R[1], Q[2], R[0], p4, p2, p4, p2, f21, &p_f21);
//  call poly_mul1(Q(:,4), f21, p4, p_f21, f22, p_f22)
  poly_mul1(Q[4], f21, p4, p_f21, f22, &p_f22);
//  call poly_sub1(f20, f22, p_f20, p_f22, f23, p_f23)
  poly_sub1(f20, f22, p_f20, p_f22, f23, &p_f23);
//  call poly_mul1(R(:,0), f23, p2, p_f23, f24, p_f24)
  poly_mul1(R[0], f23, p2, p_f23, f24, &p_f24);
//  call poly_sub1(f17, f24, p_f17, p_f24, f25, p_f25)
  poly_sub1(f17, f24, p_f17, p_f24, f25, &p_f25);
//  call poly_mul_sub1(Q(:,4), f12, R(:,2), f25, p4, p_f12, p2, p_f25, f26, p_f26)
  poly_mul_sub1(Q[4], f12, R[2], f25, p4, p_f12, p2, p_f25, f26, &p_f26);
//  call poly_mul_sub1(Q(:,0), f9, R(:,0), f26, p4, p_f9, p2, p_f26, poly_coeff, p_final)
  poly_mul_sub1(Q[0], f9, R[0], f26, p4, p_f9, p2, p_f26, poly_coeff, &p_final);

//  if (p_final /= 16) then
//     print*, 'Error. Degree of polynomial is not 16!'
//     stop
//  end if
  if (p_final != 16)
   {
    printf("Error. Degree of polynomial is not 16!\n");
    exit(1);
   }

//  if (poly_coeff(16) < 0.0d0) then
//     poly_coeff(0:16) = -poly_coeff(0:16) 
//  end if
  if (poly_coeff[16] < 0.0e0)
    for(i=0;i<17;i++)
     poly_coeff[i] *= -1.0; 

//  if (print_level > 0) then
//     print*, 'poly_coeff'
//     do i = 0, 16
//        write(*,"(i5,e15.6)") i, poly_coeff(i)
//     end do
//  end if
  if (print_level > 0)
   {
     printf("poly_coeff\n");
     for(i=0;i<17;i++)
        printf("%5d%15.6f\n", i, poly_coeff[i]);
   }

  return;
//end subroutine get_poly_coeff
 }
//!----------------------------------------------------------------------------
//subroutine poly_mul_sub2(u1, u2, u3, u4, p1, p2, p3, p4, u5, p5)
void poly_mul_sub2(double u1[5][5], double u2[5][5], double u3[5][5], double u4[5][5], int p1[2], int p2[2], int p3[2], int p4[2], double u5[5][5], int p5[2])
 {
//  implicit none
//  real(dp), dimension(0:4,0:4), intent(in) :: u1, u2, u3, u4
//  integer, dimension(2), intent(in) :: p1, p2, p3, p4
//  real(dp), dimension(0:4,0:4), intent(out) :: u5
//  integer, dimension(2), intent(out) :: p5
//  real(dp), dimension(0:4,0:4) :: d1, d2
  double d1[5][5], d2[5][5];
//  integer, dimension(2) :: pd1, pd2
  int pd1[2], pd2[2];

//  call poly_mul2(u1, u2, p1, p2, d1, pd1)
  poly_mul2(u1, u2, p1, p2, d1, pd1);
//  call poly_mul2(u3, u4, p3, p4, d2, pd2)
  poly_mul2(u3, u4, p3, p4, d2, pd2);
//  call poly_sub2(d1, d2, pd1, pd2, u5, p5)
  poly_sub2(d1, d2, pd1, pd2, u5, p5);

//end subroutine poly_mul_sub2
 }
//!----------------------------------------------------------------------------
//subroutine poly_mul2(u1, u2, p1, p2, u3, p3)
void poly_mul2(double u1[5][5], double u2[5][5], int p1[2], int p2[2], double u3[5][5], int p3[2])
 {
//  implicit none
//  real(dp), dimension(0:4,0:4), intent(in) :: u1, u2
//  integer, dimension(2), intent(in) :: p1, p2
//  real(dp), dimension(0:4,0:4), intent(out) :: u3
//  integer, intent(out) :: p3(2)
//  integer :: i1, j1, i2, j2, i3, j3, p11, p12, p21, p22
  int i1, j1, i2, j2, i3, j3, p11, p12, p21, p22;
  int i, j;
//  real(dp) :: u1ij
  double u1ij;

//  p3(:) = p1(:) + p2(:)
  for(i=0;i<2;i++)
    p3[i] = p1[i] + p2[i];
  for(i=0;i<5;i++)
    for(j=0;j<5;j++)
      u3[i][j] = 0.0e0;

//  p11 = p1(1)
  p11 = p1[0];
//  p12 = p1(2)
  p12 = p1[1];
//  p21 = p2(1)
  p21 = p2[0];
//  p22 = p2(2)
  p22 = p2[1];

//  do i1 = 0, p12
  for(i1=0;i1<=p12;i1++)
   {
//     do j1 = 0, p11
     for(j1=0;j1<=p11;j1++)
      {
//        u1ij = u1(j1,i1)
        u1ij = u1[i1][j1];
//        do i2 = 0, p22
	for(i2=0;i2<=p22;i2++)
	 {
//           i3 = i1 + i2
           i3 = i1 + i2;
//           do j2 = 0, p21
	   for (j2=0;j2<=p21;j2++)
	    {
//              j3 = j1 + j2
              j3 = j1 + j2;
//              u3(j3,i3) = u3(j3,i3) + u1ij*u2(j2,i2)
              u3[i3][j3] = u3[i3][j3] + u1ij*u2[i2][j2];
//           end do
	    }
//        end do
	 }
//     end do
      }
//  end do
   }

//end subroutine poly_mul2
 }
//!----------------------------------------------------------------------------
//subroutine poly_sub2(u1, u2, p1, p2, u3, p3)
void poly_sub2(double u1[5][5], double u2[5][5], int p1[2], int p2[2], double u3[5][5], int p3[2])
 {
//  implicit none
//  real(dp), dimension(0:4,0:4), intent(in) :: u1, u2
//  integer, intent(in) :: p1(2), p2(2)
//  real(dp), dimension(0:4,0:4), intent(out) :: u3
//  integer, intent(out) :: p3(2)
//  integer :: i, j, p11, p12, p21, p22
  int i, j, p11, p12, p21, p22;
//  logical :: i1_ok, i2_ok
  int i1_ok, i2_ok;

//  p11 = p1(1)
  p11 = p1[0];
//  p12 = p1(2)
  p12 = p1[1];
//  p21 = p2(1)
  p21 = p2[0];
//  p22 = p2(2)
  p22 = p2[1];
//  p3(1) = max(p11,p21)
  p3[0] = max(p11,p21);
//  p3(2) = max(p12,p22)
  p3[1] = max(p12,p22);

//  do i = 0, p3(2)
  for(i=0;i<=p3[1];i++)
   {
//     i1_ok = (i > p12)
    i1_ok = (i > p12);
//     i2_ok = (i > p22)
    i2_ok = (i > p22);
//     do j = 0, p3(1)
    for(j=0;j<=p3[0];j++)
     {
//        if (i2_ok .or. (j > p21)) then
//           u3(j,i) = u1(j,i)
      if (i2_ok || (j > p21))
        u3[i][j] = u1[i][j];
//        else if (i1_ok .or. (j > p11)) then
//           u3(j,i) = -u2(j,i)
      else if (i1_ok || (j > p11))
        u3[i][j] = -u2[i][j];
//        else
//           u3(j,i) = u1(j,i) - u2(j,i)
      else
        u3[i][j] = u1[i][j] - u2[i][j];
//        end if
//     end do
     }
//  end do
   }

  return;
//end subroutine poly_sub2
 }
//!----------------------------------------------------------------------------
//subroutine poly_mul_sub1(u1, u2, u3, u4, p1, p2, p3, p4, u5, p5)
void poly_mul_sub1(double u1[17], double u2[17], double u3[17], double u4[17], int p1, int p2, int p3, int p4, double u5[17], int *p5)
 {
//  implicit none
//  real(dp), dimension(0:16), intent(in) :: u1, u2, u3, u4
//  integer, intent(in) :: p1, p2, p3, p4
//  real(dp), dimension(0:16), intent(out) :: u5
//  integer, intent(out) :: p5
//  real(dp), dimension(0:16) :: d1, d2
  double d1[17], d2[17];
//  integer :: pd1, pd2
  int pd1, pd2;

//  call poly_mul1(u1, u2, p1, p2, d1, pd1)
  poly_mul1(u1, u2, p1, p2, d1, &pd1);
//  call poly_mul1(u3, u4, p3, p4, d2, pd2)
  poly_mul1(u3, u4, p3, p4, d2, &pd2);
//  call poly_sub1(d1, d2, pd1, pd2, u5, p5)
  poly_sub1(d1, d2, pd1, pd2, u5, p5);

  return;
//end subroutine poly_mul_sub1
 }
//!----------------------------------------------------------------------------
//subroutine poly_mul1(u1, u2, p1, p2, u3, p3)
void poly_mul1(double u1[17], double u2[17], int p1, int p2, double u3[17], int *p3)
 {
//  implicit none
//  real(dp), dimension(0:16), intent(in) :: u1, u2
//  integer, intent(in) :: p1, p2
//  real(dp), dimension(0:16), intent(out) :: u3
//  integer, intent(out) :: p3
//  integer :: i1, i2, i3
  int i, i1, i2, i3;
//  real(dp) :: u1i
  double u1i;

//  p3 = p1 + p2
  *p3 = p1 + p2;
//  u3(:) = 0.0d0
  for(i=0;i<17;i++)
    u3[i] = 0.;

//  do i1 = 0, p1
  for(i1=0;i1<=p1;i1++)
   {
//     u1i = u1(i1)
    u1i = u1[i1];
//     do i2 = 0, p2
    for(i2=0;i2<=p2;i2++)
     {
//        i3 = i1 + i2
        i3 = i1 + i2;
//        u3(i3) = u3(i3) + u1i*u2(i2)
        u3[i3] = u3[i3] + u1i*u2[i2];
//     end do
     }
//  end do
   }

  return;
//end subroutine poly_mul1
 }
//!----------------------------------------------------------------------------
//subroutine poly_sub1(u1, u2, p1, p2, u3, p3)
void poly_sub1(double u1[17], double u2[17], int p1, int p2, double u3[17], int *p3)
 {
//  implicit none
//  real(dp), dimension(0:16), intent(in) :: u1, u2
//  integer, intent(in) :: p1, p2
//  real(dp), dimension(0:16), intent(out) :: u3
//  integer, intent(out) :: p3
//  integer :: i
  int i;

//  p3 = max(p1, p2)
  *p3 = max(p1, p2);

//  do i = 0, p3
  for(i=0;i<=*p3;i++)
   {
//     if (i > p2) then
//        u3(i) = u1(i)
    if (i > p2)
      u3[i] = u1[i];
//     else if (i > p1) then
//        u3(i) = -u2(i)
    else if (i > p1)
      u3[i] = -u2[i];
//     else
//        u3(i) = u1(i) - u2(i)
    else
      u3[i] = u1[i] - u2[i];
//     end if
//  end do
   }

  return;
//end subroutine poly_sub1
 }
//!----------------------------------------------------------------------------
//subroutine coord_from_poly_roots(n_soln, roots, r_n1, r_a1, r_a3, r_c3, r_soln_n, r_soln_a, r_soln_c)
void coord_from_poly_roots(int *n_soln, double roots[max_soln], double r_n1[3], double r_a1[3], double r_a3[3], double r_c3[3], double r_soln_n[max_soln][3][3], double r_soln_a[max_soln][3][3], double r_soln_c[max_soln][3][3])
 {
//  implicit none
//  integer, intent(in) :: n_soln
//  real(dp), intent(in) :: r_n1(:), r_a1(:), r_a3(:), r_c3(:), roots(n_soln)
//  real(dp), intent(out) :: r_soln_n(:,:,:), r_soln_a(:,:,:), r_soln_c(:,:,:)
//  real(dp) :: ex(3), ey(3), ez(3), b_a1a2(3), b_a3a2(3), r_tmp(3)
  double ex[3], ey[3], ez[3], b_a1a2[3], b_a3a2[3], r_tmp[3];
//  real(dp) :: p_s(3,3), s1(3,3), s2(3,3), p_t(3,3), t1(3,3), t2(3,3)
  double p_s[3][3], s1[3][3], s2[3][3], p_t[3][3], t1[3][3], t2[3][3];
//  real(dp) :: p_s_c(3,3), s1_s(3,3), s2_s(3,3), p_t_c(3,3), t1_s(3,3), t2_s(3,3)
  double p_s_c[3][3], s1_s[3][3], s2_s[3][3], p_t_c[3][3], t1_s[3][3], t2_s[3][3];
//  real(dp) :: angle, sig1_init, half_tan(3)
  double angle, sig1_init, half_tan[3];
//  real(dp) :: cos_tau(0:3), sin_tau(0:3), cos_sig(3), sin_sig(3), ht, tmp, sig1
  double cos_tau[4], sin_tau[4], cos_sig[3], sin_sig[3], ht, tmp, sig1;
//  real(dp) :: r_s(3), r_t(3), r0(3), r_n(3,3), r_a(3,3), r_c(3,3), p(4), Us(3,3)
  double r_s[3], r_t[3], r0[3], r_n[3][3], r_a[3][3], r_c[3][3], p[4], Us[3][3];
//  integer :: i_soln, i, j
  int i_soln, i, j;
//  real(dp) :: a1c1, c1n2, n2a2, a2c2, c2n3, n3a3, a1a2, a2a3
  double a1c1, c1n2, n2a2, a2c2, c2n3, n3a3, a1a2, a2a3;
//  real(dp) :: rr_a1c1(3), rr_c1n2(3), rr_n2a2(3), rr_a2c2(3), rr_c2n3(3), rr_n3a3(3), rr_a1a2(3), rr_a2a3(3)
  double rr_a1c1[3], rr_c1n2[3], rr_n2a2[3], rr_a2c2[3], rr_c2n3[3], rr_n3a3[3], rr_a1a2[3], rr_a2a3[3];
//  real(dp) :: a3a1a2, a2a3a1, n1a1c1, n2a2c2, n3a3c3, a1c1n2a2, a2c2n3a3
  double a3a1a2, a2a3a1, n1a1c1, n2a2c2, n3a3c3, a1c1n2a2, a2c2n3a3;
  double tmp_value, ex_tmp[3];
  double tmp_array[3], tmp_array1[3], tmp_array2[3], tmp_array3[3];
  double mat1[3], mat2[3], mat3[3], mat4[3], mat5[3];
  double mat11[3], mat22[3], mat33[3], mat44[3], mat55[3];
  
//  if (n_soln == 0) return
  if (*n_soln == 0)
    return;

//  ! Define body frame (ex, ey, ez)
//  ex(:) = b_a1a3(:)
  for(i=0;i<3;i++)
    ex[i] = b_a1a3[i];
//  call cross(r_a1n1, ex, ez)
  cross(r_a1n1, ex, ez);
//  ez(:) = ez(:)/sqrt(dot_product(ez,ez))
  tmp_value = sqrt(dot_product(ez,ez));
  for(i=0;i<3;i++)
    ez[i] = ez[i]/tmp_value;
//  call cross(ez, ex, ey)
  cross(ez, ex, ey);
//  ! vertual bond vectors in the reference plane
//  b_a1a2(:) = -cos_alpha(1)*ex(:) + sin_alpha(1)*ey(:)
//  b_a3a2(:) = cos_alpha(3)*ex(:) + sin_alpha(3)*ey(:)
  for(i=0;i<3;i++)
   {
    b_a1a2[i] = -cos_alpha[0]*ex[i] + sin_alpha[0]*ey[i];
    b_a3a2[i] = cos_alpha[2]*ex[i] + sin_alpha[2]*ey[i];
   }
//  !! Define cone coordinates for each angle joint.
//  ! (p_s,s1,s2) and (p_t,t1,t2):  Right Orthonormal systems
//  ! residue 1
//  p_s(:,1) = -ex(:)
//  s1(:,1)  = ez(:)  ! (p_s)X(p_t)/||(p_s)X(p_t)||
//  s2(:,1)  = ey(:)  ! p_s X s1
//  p_t(:,1) = b_a1a2(:)
//  t1(:,1)  = ez(:)  ! s1
//  t2(:,1)  = sin_alpha(1)*ex(:) + cos_alpha(1)*ey(:) ! p_t X t1
  for(i=0;i<3;i++)
   {
    p_s[0][i] = -ex[i];
    s1[0][i]  = ez[i];
    s2[0][i]  = ey[i];
    p_t[0][i] = b_a1a2[i];
    t1[0][i]  = ez[i];
    t2[0][i]  = sin_alpha[0]*ex[i] + cos_alpha[0]*ey[i];
   }
//  ! residue 2
//  p_s(:,2) = -b_a1a2(:)
//  s1(:,2)  = -ez(:)
//  s2(:,2)  = t2(:,1)  ! sina1*ex(:) + cosa1*ey(:)
//  p_t(:,2) = -b_a3a2(:)
//  t1(:,2)  = -ez(:)
//  t2(:,2)  = sin_alpha(3)*ex(:) - cos_alpha(3)*ey(:) 
  for(i=0;i<3;i++)
   {
    p_s[1][i] = -b_a1a2[i];
    s1[1][i]  = -ez[i];
    s2[1][i]  = t2[0][i];
    p_t[1][i] = -b_a3a2[i];
    t1[1][i]  = -ez[i];
    t2[1][i]  = sin_alpha[2]*ex[i] - cos_alpha[2]*ey[i];
   }
//  ! residue 3
//  p_s(:,3) = b_a3a2(:)
//  s2(:,3)  = t2(:,2)   ! sina3*ex(:) + cosa3*ey(:) 
//  s1(:,3)  = ez(:)  
//  p_t(:,3) = ex(:)
//  t1(:,3) =  ez(:) 
//  t2(:,3) = -ey(:) 
  for(i=0;i<3;i++)
   {
    p_s[2][i] = b_a3a2[i];
    s2[2][i]  = t2[1][i]; 
    s1[2][i]  = ez[i];
    p_t[2][i] = ex[i];
    t1[2][i] =  ez[i];
    t2[2][i] = -ey[i];
   }
//  ! scale vectors
//  do i = 1, 3
//     p_s_c(:,i) = p_s(:,i)*cos_xi(i)
//     s1_s(:,i)  =  s1(:,i)*sin_xi(i)
//     s2_s(:,i)  =  s2(:,i)*sin_xi(i)
//     p_t_c(:,i) = p_t(:,i)*cos_eta(i)
//     t1_s(:,i)  =  t1(:,i)*sin_eta(i)
//     t2_s(:,i)  =  t2(:,i)*sin_eta(i)
//  end do
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
     {
      p_s_c[i][j] = p_s[i][j]*cos_xi[i];
      s1_s[i][j]  = s1[i][j]*sin_xi[i];
      s2_s[i][j]  = s2[i][j]*sin_xi[i];
      p_t_c[i][j] = p_t[i][j]*cos_eta[i];
      t1_s[i][j]  = t1[i][j]*sin_eta[i];
      t2_s[i][j]  = t2[i][j]*sin_eta[i];
     }

//  ! initial sig(1)
//  r_tmp(:) = (r_a1n1(:)/len_na(1) - p_s_c(:,1))/sin_xi(1)
  for(i=0;i<3;i++)
    r_tmp[i] = (r_a1n1[i]/len_na[0] - p_s_c[0][i])/sin_xi[0];
//  call calc_bnd_ang(s1(:,1), r_tmp, angle)
  calc_bnd_ang(s1[0], r_tmp, &angle);
//  sig1_init = sign(angle, dot_product(r_tmp(:),s2(:,1)))
  sig1_init = sign(angle, dot_product(r_tmp,s2[0]));

//  ! CA
//  r_a(:,1) = r_a1(:)
//  r_a(:,2) = r_a1(:) + len_aa(2)*b_a1a2(:)
//  r_a(:,3) = r_a3(:)
//  r0(:) = r_a1(:)
  for(i=0;i<3;i++)
   {
    r_a[0][i] = r_a1[i];
    r_a[1][i] = r_a1[i] + len_aa[1]*b_a1a2[i];
    r_a[2][i] = r_a3[i];
    r0[i] = r_a1[i];
   }

//  do i_soln = 1, n_soln 
  for(i_soln=0;i_soln<*n_soln;i_soln++)
   {
//     half_tan(3) = roots(i_soln)
     half_tan[2] = roots[i_soln];
//     half_tan(2) = calc_t2(half_tan(3))
     half_tan[1] = calc_t2(half_tan[2]);
//     half_tan(1) = calc_t1(half_tan(3), half_tan(2))
     half_tan[0] = calc_t1(half_tan[2], half_tan[1]);
//     do i = 1, 3
//        ht = half_tan(i)
//        tmp = 1.0d0 + ht*ht
//        cos_tau(i) = (1.0d0 - ht*ht)/tmp
//        sin_tau(i) = 2.0d0*ht/tmp
//     end do
     for(i=1;i<=3;i++)
      {
       ht = half_tan[i-1];
       tmp = 1.0e0 + ht*ht;
       cos_tau[i] = (1.0e0 - ht*ht)/tmp;
       sin_tau[i] = 2.0e0*ht/tmp;
      }
//     cos_tau(0) = cos_tau(3)
     cos_tau[0] = cos_tau[3];
//     sin_tau(0) = sin_tau(3)
     sin_tau[0] = sin_tau[3];
//     printf("cos: %7.3f%7.3f%7.3f%7.3f\n", cos_tau[0], cos_tau[1], cos_tau[2], cos_tau[3]);
//     printf("sin: %7.3f%7.3f%7.3f%7.3f\n", sin_tau[0], sin_tau[1], sin_tau[2], sin_tau[3]);
//     do i = 1, 3
//        j = i - 1
//        cos_sig(i) = cos_delta(j)*cos_tau(j) + sin_delta(j)*sin_tau(j)
//        sin_sig(i) = sin_delta(j)*cos_tau(j) - cos_delta(j)*sin_tau(j)
//     end do
     for(i=0;i<3;i++)
      {
       cos_sig[i] = cos_delta[i]*cos_tau[i] + sin_delta[i]*sin_tau[i];
       sin_sig[i] = sin_delta[i]*cos_tau[i] - cos_delta[i]*sin_tau[i];
      }
//     do i = 1, 3
//        r_s(:) = p_s_c(:,i) + cos_sig(i)*s1_s(:,i) + sin_sig(i)*s2_s(:,i)
//        r_t(:) = p_t_c(:,i) + cos_tau(i)*t1_s(:,i) + sin_tau(i)*t2_s(:,i) 
//        r_n(:,i) = r_s(:)*len_na(i) + r_a(:,i)
//        r_c(:,i) = r_t(:)*len_ac(i) + r_a(:,i)
//     end do
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
       {
        r_s[j] = p_s_c[i][j] + cos_sig[i]*s1_s[i][j] + sin_sig[i]*s2_s[i][j];
        r_t[j] = p_t_c[i][j] + cos_tau[i+1]*t1_s[i][j] + sin_tau[i+1]*t2_s[i][j];
        r_n[i][j] = r_s[j]*len_na[i] + r_a[i][j];
        r_c[i][j] = r_t[j]*len_ac[i] + r_a[i][j];
       }
      
//     ! rotate back atoms by -(sig(1) - sig1_init) around -ex
//     sig1 = atan2(sin_sig(1), cos_sig(1))
     sig1 = atan2(sin_sig[0], cos_sig[0]);
//     call quaternion(-ex, -(sig1 - sig1_init)*0.25d0, p)
     ex_tmp[0] = -ex[0];
     ex_tmp[1] = -ex[1];
     ex_tmp[2] = -ex[2];
     tmp_value = -(sig1-sig1_init)*0.25;
     quaternion(ex_tmp, tmp_value, p);
//     call rotation_matrix(p, Us)
     rotation_matrix(p, Us);
     
//     r_soln_n(:,1,i_soln) = r_n1(:)
//     r_soln_a(:,1,i_soln) = r_a1(:)
//     r_soln_c(:,1,i_soln) = matmul(Us, r_c(:,1) - r0(:)) + r0(:)
//     r_soln_n(:,2,i_soln) = matmul(Us, r_n(:,2) - r0(:)) + r0(:)
//     r_soln_a(:,2,i_soln) = matmul(Us, r_a(:,2) - r0(:)) + r0(:)
//     r_soln_c(:,2,i_soln) = matmul(Us, r_c(:,2) - r0(:)) + r0(:)
//     r_soln_n(:,3,i_soln) = matmul(Us, r_n(:,3) - r0(:)) + r0(:)
//     r_soln_a(:,3,i_soln) = r_a3(:)
//     r_soln_c(:,3,i_soln) = r_c3(:)
     for(i=0;i<3;i++)
      {
       mat11[i] = r_c[0][i]-r0[i];
       mat22[i] = r_n[1][i]-r0[i];
       mat33[i] = r_a[1][i]-r0[i];
       mat44[i] = r_c[1][i]-r0[i];
       mat55[i] = r_n[2][i]-r0[i];
      }
     matmul(Us,mat11,mat1);
     matmul(Us,mat22,mat2);
     matmul(Us,mat33,mat3);
     matmul(Us,mat44,mat4);
     matmul(Us,mat55,mat5);
     for(i=0;i<3;i++)
      {
       r_soln_n[i_soln][0][i] = r_n1[i];
       r_soln_a[i_soln][0][i] = r_a1[i];
       r_soln_c[i_soln][0][i] = mat1[i] + r0[i];
       r_soln_n[i_soln][1][i] = mat2[i] + r0[i];
       r_soln_a[i_soln][1][i] = mat3[i] + r0[i];
       r_soln_c[i_soln][1][i] = mat4[i] + r0[i];
       r_soln_n[i_soln][2][i] = mat5[i] + r0[i];
       r_soln_a[i_soln][2][i] = r_a3[i];
       r_soln_c[i_soln][2][i] = r_c3[i];
      }


//     if (print_level > 0) then
     if (print_level > 0)
      {
//        print*, 'roots: t0, t2, t1', i_soln
        printf("roots: t0, t2, t1 %d\n", i_soln);
//        write(*,"(3f15.6)") half_tan(3), half_tan(2), half_tan(1)
        printf("%15.6f %15.6f %15.6f\n", half_tan[2], half_tan[1], half_tan[0]);

//        rr_a1c1(:) = r_soln_c(:,1,i_soln) - r_soln_a(:,1,i_soln)
//        rr_c1n2(:) = r_soln_n(:,2,i_soln) - r_soln_c(:,1,i_soln)
//        rr_n2a2(:) = r_soln_a(:,2,i_soln) - r_soln_n(:,2,i_soln)
//        rr_a2c2(:) = r_soln_c(:,2,i_soln) - r_soln_a(:,2,i_soln)
//        rr_c2n3(:) = r_soln_n(:,3,i_soln) - r_soln_c(:,2,i_soln)
//        rr_n3a3(:) = r_soln_a(:,3,i_soln) - r_soln_n(:,3,i_soln)
//        rr_a1a2(:) = r_soln_a(:,2,i_soln) - r_soln_a(:,1,i_soln)
//        rr_a2a3(:) = r_soln_a(:,3,i_soln) - r_soln_a(:,2,i_soln)
	for(i=0;i<3;i++)
	 {
          rr_a1c1[i] = r_soln_c[i_soln][0][i] - r_soln_a[i_soln][0][i];
          rr_c1n2[i] = r_soln_n[i_soln][1][i] - r_soln_c[i_soln][0][i];
          rr_n2a2[i] = r_soln_a[i_soln][1][i] - r_soln_n[i_soln][1][i];
          rr_a2c2[i] = r_soln_c[i_soln][1][i] - r_soln_a[i_soln][1][i];
          rr_c2n3[i] = r_soln_n[i_soln][2][i] - r_soln_c[i_soln][1][i];
          rr_n3a3[i] = r_soln_a[i_soln][2][i] - r_soln_n[i_soln][2][i];
          rr_a1a2[i] = r_soln_a[i_soln][1][i] - r_soln_a[i_soln][0][i];
          rr_a2a3[i] = r_soln_a[i_soln][2][i] - r_soln_a[i_soln][1][i];
	 }
	      
//        a1c1 = sqrt(dot_product(rr_a1c1, rr_a1c1))
//        c1n2 = sqrt(dot_product(rr_c1n2, rr_c1n2))
//        n2a2 = sqrt(dot_product(rr_n2a2, rr_n2a2))
//        a2c2 = sqrt(dot_product(rr_a2c2, rr_a2c2))
//        c2n3 = sqrt(dot_product(rr_c2n3, rr_c2n3))
//        n3a3 = sqrt(dot_product(rr_n3a3, rr_n3a3))
//       a1a2 = sqrt(dot_product(rr_a1a2, rr_a1a2))
//        a2a3 = sqrt(dot_product(rr_a2a3, rr_a2a3))
        a1c1 = sqrt(dot_product(rr_a1c1, rr_a1c1));
        c1n2 = sqrt(dot_product(rr_c1n2, rr_c1n2));
        n2a2 = sqrt(dot_product(rr_n2a2, rr_n2a2));
        a2c2 = sqrt(dot_product(rr_a2c2, rr_a2c2));
        c2n3 = sqrt(dot_product(rr_c2n3, rr_c2n3));
        n3a3 = sqrt(dot_product(rr_n3a3, rr_n3a3));
        a1a2 = sqrt(dot_product(rr_a1a2, rr_a1a2));
        a2a3 = sqrt(dot_product(rr_a2a3, rr_a2a3));

//        write(*,"('na: n2a2, n3a3 = ', 4f9.3)") len0(3), n2a2, len0(6), n3a3
        printf("na: n2a2, n3a3 = %9.3f%9.3f%9.3f%9.3f\n", len0[2], n2a2, len0[5], n3a3);
//        write(*,"('ac: a1c1, a2c2 = ', 4f9.3)") len0(1), a1c1, len0(4), a2c2
        printf("ac: a1c1, a2c2 = %9.3f%9.3f%9.3f%9.3f\n", len0[0], a1c1, len0[3], a2c2);
//        write(*,"('cn: c1n2, c2n3 = ', 4f9.3)") len0(2), c1n2, len0(5), c2n3
        printf("cn: c1n2, c2n3 = %9.3f%9.3f%9.3f%9.3f\n", len0[1], c1n2, len0[4], c2n3);
//        write(*,"('aa: a1a2, a2a3 = ', 4f9.3)") len_aa(2), a1a2, len_aa(3), a2a3
        printf("aa: a1a2, a2a3 = %9.3f%9.3f%9.3f%9.3f\n", len_aa[1], a1a2, len_aa[2], a2a3);

//        call calc_bnd_ang(b_a1a3, rr_a1a2/a1a2, a3a1a2)
	for(i=0;i<3;i++)
	  tmp_array[i] = rr_a1a2[i]/a1a2;
        calc_bnd_ang(b_a1a3, tmp_array, &a3a1a2);
//        call calc_bnd_ang(rr_a2a3/a2a3, b_a1a3, a2a3a1)
	for(i=0;i<3;i++)
	  tmp_array[i] = rr_a2a3[i]/a2a3;
        calc_bnd_ang(tmp_array, b_a1a3, &a2a3a1);
//        write(*,"('alpha1, alpha3 = ', 2f9.3)") (pi-a3a1a2)*rad2deg, (pi-a2a3a1)*rad2deg
        printf("alpha1, alpha3 = %9.3f%9.3f\n", (pi-a3a1a2)*rad2deg, (pi-a2a3a1)*rad2deg);
//        call calc_bnd_ang(b_a1n1, rr_a1c1/a1c1, n1a1c1)
	for(i=0;i<3;i++)
	  tmp_array[i] = rr_a1c1[i]/a1c1;
        calc_bnd_ang(b_a1n1, tmp_array, &n1a1c1);
//        call calc_bnd_ang(b_a3c3, -rr_n3a3/n3a3, n3a3c3)
	for(i=0;i<3;i++)
	  tmp_array[i] = -rr_n3a3[i]/n3a3;
        calc_bnd_ang(b_a3c3, tmp_array, &n3a3c3);
//        call calc_bnd_ang(rr_a2c2/a2c2, -rr_n2a2/n2a2, n2a2c2)
	for(i=0;i<3;i++)
	 {
	  tmp_array1[i] = rr_a2c2[i]/a2c2;
	  tmp_array2[i] = -rr_n2a2[i]/n2a2;
	 }
        calc_bnd_ang(tmp_array1, tmp_array2, &n2a2c2);
//        write(*,"('ang_nac = ', 2f9.3)") b_ang0(1)*rad2deg, n1a1c1*rad2deg
        printf("ang_nac = %9.3f%9.3f\n", b_ang0[0]*rad2deg, n1a1c1*rad2deg);
//        write(*,"('ang_nac = ', 2f9.3)") b_ang0(4)*rad2deg, n2a2c2*rad2deg 
        printf("ang_nac = %9.3f%9.3f\n", b_ang0[3]*rad2deg, n2a2c2*rad2deg);
//        write(*,"('ang_nac = ', 2f9.3)") b_ang0(7)*rad2deg, n3a3c3*rad2deg
        printf("ang_nac = %9.3f%9.3f\n", b_ang0[6]*rad2deg, n3a3c3*rad2deg);

//        call calc_dih_ang(rr_a1c1/a1c1, rr_c1n2/c1n2, rr_n2a2/n2a2, a1c1n2a2)
	for(i=0;i<3;i++)
	 {
	  tmp_array1[i] = rr_a1c1[i]/a1c1;
	  tmp_array2[i] = rr_c1n2[i]/c1n2;
	  tmp_array3[i] = rr_n2a2[i]/n2a2;
	 }
        calc_dih_ang(tmp_array1, tmp_array2, tmp_array3, &a1c1n2a2);
//        call calc_dih_ang(rr_a2c2/a2c2, rr_c2n3/c2n3, rr_n3a3/n3a3, a2c2n3a3)
	for(i=0;i<3;i++)
	 {
	  tmp_array1[i] = rr_a2c2[i]/a2c2;
	  tmp_array2[i] = rr_c2n3[i]/c2n3;
	  tmp_array3[i] = rr_n3a3[i]/n3a3;
	 }
        calc_dih_ang(tmp_array1, tmp_array2, tmp_array3, &a2c2n3a3);
//        write(*,"('t_ang1 = ', 2f9.3)") t_ang0(1)*rad2deg, a1c1n2a2*rad2deg
        printf("t_ang1 = %9.3f%9.3f\n", t_ang0[0]*rad2deg, a1c1n2a2*rad2deg);
//        write(*,"('t_ang2 = ', 2f9.3)") t_ang0(2)*rad2deg, a2c2n3a3*rad2deg
        printf("t_ang2 = %9.3f%9.3f\n", t_ang0[1]*rad2deg, a2c2n3a3*rad2deg);
//     end if
      }

//  end do
   }

  return;
//end subroutine coord_from_poly_roots
 }
//!-----------------------------------------------------------------------
//function calc_t2(t0)
double calc_t2(double t0)
 {
//  implicit none
//  real(dp), intent(in) :: t0
//  real(dp) :: calc_t2
  double tmp_value;
//  real(dp) :: B0, B1, B2, A0, A1, A2, A3, A4, B2_2, B2_3
  double B0, B1, B2, A0, A1, A2, A3, A4, B2_2, B2_3;
//  real(dp) :: K0, K1, K2, K3, t0_2, t0_3, t0_4
  double K0, K1, K2, K3, t0_2, t0_3, t0_4;

//  t0_2 = t0*t0
  t0_2 = t0*t0;
//  t0_3 = t0_2*t0
  t0_3 = t0_2*t0;
//  t0_4 = t0_3*t0
  t0_4 = t0_3*t0;

//  A0 = Q(0,0) + Q(1,0)*t0 + Q(2,0)*t0_2 + Q(3,0)*t0_3 + Q(4,0)*t0_4
  A0 = Q[0][0] + Q[0][1]*t0 + Q[0][2]*t0_2 + Q[0][3]*t0_3 + Q[0][4]*t0_4;
//  A1 = Q(0,1) + Q(1,1)*t0 + Q(2,1)*t0_2 + Q(3,1)*t0_3 + Q(4,1)*t0_4
  A1 = Q[1][0] + Q[1][1]*t0 + Q[1][2]*t0_2 + Q[1][3]*t0_3 + Q[1][4]*t0_4;
//  A2 = Q(0,2) + Q(1,2)*t0 + Q(2,2)*t0_2 + Q(3,2)*t0_3 + Q(4,2)*t0_4
  A2 = Q[2][0] + Q[2][1]*t0 + Q[2][2]*t0_2 + Q[2][3]*t0_3 + Q[2][4]*t0_4;
//  A3 = Q(0,3) + Q(1,3)*t0 + Q(2,3)*t0_2 + Q(3,3)*t0_3 + Q(4,3)*t0_4
  A3 = Q[3][0] + Q[3][1]*t0 + Q[3][2]*t0_2 + Q[3][3]*t0_3 + Q[3][4]*t0_4;
//  A4 = Q(0,4) + Q(1,4)*t0 + Q(2,4)*t0_2 + Q(3,4)*t0_3 + Q(4,4)*t0_4
  A4 = Q[4][0] + Q[4][1]*t0 + Q[4][2]*t0_2 + Q[4][3]*t0_3 + Q[4][4]*t0_4;

//  B0 = R(0,0) + R(1,0)*t0 + R(2,0)*t0_2
  B0 = R[0][0] + R[0][1]*t0 + R[0][2]*t0_2;
//  B1 = R(0,1) + R(1,1)*t0 + R(2,1)*t0_2
  B1 = R[1][0] + R[1][1]*t0 + R[1][2]*t0_2;
//  B2 = R(0,2) + R(1,2)*t0 + R(2,2)*t0_2
  B2 = R[2][0] + R[2][1]*t0 + R[2][2]*t0_2;

//  B2_2 = B2*B2
  B2_2 = B2*B2;
//  B2_3 = B2_2*B2
  B2_3 = B2_2*B2;

//  K0 = A2*B2 - A4*B0
  K0 = A2*B2 - A4*B0;
//  K1 = A3*B2 - A4*B1
  K1 = A3*B2 - A4*B1;
//  K2 = A1*B2_2 - K1*B0
  K2 = A1*B2_2 - K1*B0;
//  K3 = K0*B2 - K1*B1
  K3 = K0*B2 - K1*B1;
  
//  calc_t2 = (K3*B0 - A0*B2_3)/(K2*B2 - K3*B1)
  tmp_value = (K3*B0 - A0*B2_3)/(K2*B2 - K3*B1);

  return tmp_value;
//end function calc_t2
 }
//!-----------------------------------------------------------------------
//function calc_t1(t0, t2)
double calc_t1(double t0, double t2)
 {
//  implicit none
//  real(dp), intent(in) :: t0, t2
//  real(dp) :: calc_t1
  double tmp_value;
//  real(dp) :: U11, U12, U13, U31, U32, U33
  double U11, U12, U13, U31, U32, U33;
//  real(dp) :: t0_2, t2_2
  double t0_2, t2_2;

//  t0_2 = t0*t0
  t0_2 = t0*t0;
//  t2_2 = t2*t2
  t2_2 = t2*t2;

//  U11 = C0(0,1) + C0(1,1)*t0 + C0(2,1)*t0_2
  U11 = C0[0][0] + C0[0][1]*t0 + C0[0][2]*t0_2;
//  U12 = C1(0,1) + C1(1,1)*t0 + C1(2,1)*t0_2
  U12 = C1[0][0] + C1[0][1]*t0 + C1[0][2]*t0_2;
//  U13 = C2(0,1) + C2(1,1)*t0 + C2(2,1)*t0_2
  U13 = C2[0][0] + C2[0][1]*t0 + C2[0][2]*t0_2;
//  U31 = C0(0,2) + C0(1,2)*t2 + C0(2,2)*t2_2
  U31 = C0[1][0] + C0[1][1]*t2 + C0[1][2]*t2_2;
//  U32 = C1(0,2) + C1(1,2)*t2 + C1(2,2)*t2_2
  U32 = C1[1][0] + C1[1][1]*t2 + C1[1][2]*t2_2;
//  U33 = C2(0,2) + C2(1,2)*t2 + C2(2,2)*t2_2
  U33 = C2[1][0] + C2[1][1]*t2 + C2[1][2]*t2_2;

//  calc_t1 = (U31*U13-U11*U33)/(U12*U33-U13*U32)
  tmp_value = (U31*U13-U11*U33)/(U12*U33-U13*U32);

  return tmp_value;
//end function calc_t1
 }
//!-----------------------------------------------------------------------
//subroutine calc_dih_ang(r1, r2, r3, angle)
void calc_dih_ang(double r1[3], double r2[3], double r3[3], double *angle)
 {
//!-----------------------------------------------------------------------
//! r1=Rab, r2=Rbc, r3=Rcd : angle between planes abc and bcd
//!-----------------------------------------------------------------------
//  implicit none
//  real(dp), intent(in) :: r1(3), r2(3), r3(3)
//  real(dp), intent(out) :: angle
//  real(dp), dimension(3) :: p, q, s
  double p[3], q[3], s[3];
//  real(dp) :: arg
  double arg;

//  call cross(r1, r2, p)
  cross(r1, r2, p);
//  call cross(r2, r3, q)
  cross(r2, r3, q);
//  call cross(r3, r1, s)
  cross(r3, r1, s);
//  arg = dot_product(p,q)/sqrt(dot_product(p,p)*dot_product(q,q))
  arg = dot_product(p,q)/sqrt(dot_product(p,p)*dot_product(q,q));
//  arg = sign(min(abs(arg),1.0d0),arg) ! to be sure abs(arg)<=1
  arg = sign(min(fabs(arg),1.0e0),arg);
//  angle = sign(acos(arg), dot_product(s,r2))
  *angle = sign(acos(arg), dot_product(s,r2));

  return;
// end subroutine calc_dih_ang
 }
//!-----------------------------------------------------------------------
//subroutine calc_bnd_ang(r1, r2, angle)
void calc_bnd_ang(double r1[3], double r2[3], double *angle)
 {
//!-----------------------------------------------------------------------
//! assume that each vector is normalized
//! r1=Rba, r2=Rbc: angle between r1 and r2
//!-----------------------------------------------------------------------
//  implicit none
//  real(dp), intent(in) :: r1(3), r2(3)
//  real(dp), intent(out) :: angle
//  real(dp) :: arg
  double arg;

//  arg = dot_product(r1, r2)
  arg = dot_product(r1, r2);
//  arg = sign(min(abs(arg),1.0d0),arg) ! to be sure abs(arg)<=1
  arg = sign(min(fabs(arg),1.0e0),arg);
//  angle = acos(arg)
  *angle = acos(arg);
  
  return;
//end subroutine calc_bnd_ang
 }
//!-----------------------------------------------------------------------
//subroutine cross(p, q, s)
void cross(double p[3], double q[3], double s[3])
 {
//!-----------------------------------------------------------------------
//  implicit none
//  real(dp), dimension(:), intent(in) :: p, q
//  real(dp), dimension(:), intent(out) :: s

//  s(1) = p(2)*q(3) - p(3)*q(2)
  s[0] = p[1]*q[2] - p[2]*q[1];
//  s(2) = p(3)*q(1) - p(1)*q(3)
  s[1] = p[2]*q[0] - p[0]*q[2];
//  s(3) = p(1)*q(2) - p(2)*q(1)
  s[2] = p[0]*q[1] - p[1]*q[0];

  return;
//end subroutine cross
 }
//!-----------------------------------------------------------------------
//subroutine quaternion(axis, quarter_ang, p)
void quaternion(double axis[3], double quarter_ang, double p[4])
 {
//!-----------------------------------------------------------------------
//! calculate quaternion, given rotation axis and angle. 
//!-----------------------------------------------------------------------
//  implicit none
//  real(dp), dimension(:), intent(in) :: axis
//  real(dp), intent(in) :: quarter_ang
//  real(dp), dimension(:), intent(out) :: p
//  real(dp) :: tan_w, tan_sqr, tan1, cosine, sine
  double tan_w, tan_sqr, tan1, cosine, sine;

//  tan_w = tan(quarter_ang)
  tan_w = tan(quarter_ang);
//  tan_sqr = tan_w * tan_w
  tan_sqr = tan_w * tan_w;
//  tan1 = 1.0d0 + tan_sqr
  tan1 = 1.0e0 + tan_sqr;
//  cosine = (1.0d0 - tan_sqr)/tan1
  cosine = (1.0e0 - tan_sqr)/tan1;
//  sine = 2.0d0*tan_w/tan1
  sine = 2.0e0*tan_w/tan1;
//  p(1) = cosine
  p[0] = cosine;
//  p(2:4) = axis(1:3) * sine
  p[1] = axis[0] * sine;
  p[2] = axis[1] * sine;
  p[3] = axis[2] * sine;

  return;
//end subroutine quaternion
 }
//!-----------------------------------------------------------------------
//subroutine rotation_matrix(q, U)
void rotation_matrix(double q[4], double U[3][3])
 {
//!-----------------------------------------------------------------------
//! constructs rotation matrix U from quaternion q.
//!-----------------------------------------------------------------------
//  implicit none
//  real(dp), dimension(:), intent(in) :: q
//  real(dp), dimension(:,:), intent(out) :: U
//  real(dp) :: q0,q1,q2,q3,b0,b1,b2,b3,q00,q01,q02,q03,q11,q12,q13,q22,q23,q33  
  double q0,q1,q2,q3,b0,b1,b2,b3,q00,q01,q02,q03,q11,q12,q13,q22,q23,q33;

//  q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4)
  q0 = q[0];
  q1 = q[1];
  q2 = q[2];
  q3 = q[3];
//  b0 = 2.0d0*q0; b1 = 2.0d0*q1
  b0 = 2.0e0*q0;
  b1 = 2.0e0*q1;
//  q00 = b0*q0-1.0d0; q02 = b0*q2; q03 = b0*q3
  q00 = b0*q0-1.0e0;
  q02 = b0*q2;
  q03 = b0*q3;
//  q11 = b1*q1;       q12 = b1*q2; q13 = b1*q3  
  q11 = b1*q1;
  q12 = b1*q2;
  q13 = b1*q3;
//  b2 = 2.0d0*q2; b3 = 2.0d0*q3
  b2 = 2.0e0*q2;
  b3 = 2.0e0*q3;
//  q01 = b0*q1; q22 = b2*q2; q23 = b2*q3; q33 = b3*q3 
  q01 = b0*q1;
  q22 = b2*q2;
  q23 = b2*q3;
  q33 = b3*q3;
//  U(1,1) = q00+q11; U(1,2) = q12-q03; U(1,3) = q13+q02
  U[0][0] = q00+q11;
  U[0][1] = q12-q03;
  U[0][2] = q13+q02;
//  U(2,1) = q12+q03; U(2,2) = q00+q22; U(2,3) = q23-q01
  U[1][0] = q12+q03;
  U[1][1] = q00+q22;
  U[1][2] = q23-q01;
//  U(3,1) = q13-q02; U(3,2) = q23+q01; U(3,3) = q00+q33
  U[2][0] = q13-q02;
  U[2][1] = q23+q01;
  U[2][2] = q00+q33;

  return;
//end subroutine rotation_matrix
 }
//!----------------------------------------------------------------------------
//END MODULE tripep_closure
//!----------------------------------------------------------------------------
}
