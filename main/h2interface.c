/*
	tetmag - A general-purpose finite-element micromagnetic simulation software package
	Copyright (C) 2016-2023 CNRS and Université de Strasbourg

	Author: Riccardo Hertel

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU Affero General Public License as
	published by the Free Software Foundation, either version 3 of the
	License, or (at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU Affero General Public License for more details.

	Contact: Riccardo Hertel, IPCMS Strasbourg, 23 rue du Loess,
			 67034 Strasbourg, France.
		 riccardo.hertel@ipcms.unistra.fr

	You should have received a copy of the GNU Affero General Public License
	along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

/*
 * h2interface.c
 *
 *  Created on: Jul 2, 2017
 *      Author: riccardo
 */

#include "basic.h"
#include "laplacebem3d.h"
#include "string.h"
#include "matrixnorms.h"
#include "LindholmCWrapper.h"

static ph2matrix Kh2;
static pavector b;
static pavector diagPart;
// static phmatrix HL;
static psurface3d gr2;
int readingH2;
static const char *filenameC;

int getNumberOfVertices()
{
	return gr2->vertices;
}

void readH2MatrixFromFile(int read, const char *filename)
{
	readingH2 = read;
	filenameC = filename;
}

void deleteGeometry()
{
	del_surface3d(gr2);
}

void defineGeometry(int bnx, int nbel, int ned, double *xyz, int *bij, double *surf, double *nv, int *el_ed, int *ed)
{

	gr2 = new_surface3d(bnx, ned, nbel);
	gr2->triangles = nbel;
	gr2->vertices = bnx;

	memcpy(gr2->x, xyz, 3 * bnx * sizeof(double));
	memcpy(gr2->t, bij, 3 * nbel * sizeof(uint));
	memcpy(gr2->g, surf, nbel * sizeof(double));
	memcpy(gr2->n, nv, 3 * nbel * sizeof(double));
	memcpy(gr2->s, el_ed, 3 * nbel * sizeof(uint));
	memcpy(gr2->e, ed, 2 * ned * sizeof(uint));
}

double *H2_mvp_sub(double *v, int dim)
{
	pavector f;
	f = new_pointer_avector(v, dim);
	clear_avector(b);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (uint i = 0; i < b->dim; ++i)
	{
		b->v[i] = f->v[i] * diagPart->v[i];
	}

	addeval_h2matrix_avector(1.0, (ph2matrix)Kh2, f, b);
	// ^^^ This calculates b <- 1.0 * Kh2 * f + b

	del_avector(f);
	return b->v;
}

double *H2_mvp(double *v, int dim)
{
	pavector f;
	f = new_pointer_avector(v, dim);
	clear_avector(b);

	addeval_h2matrix_avector(1.0, (ph2matrix)Kh2, f, b);

	del_avector(f);
	return b->v;
}

/*
double* HLibMVP(double* v, int dim) {
	pavector f;
	f = new_pointer_avector(v, dim);

	clear_avector(b);
	addeval_hmatrix_avector(1.0, (phmatrix) HL, f, b);

	del_avector(f);
	return b->v;
}
*/

// ##################################### begin collocation ##################################################"

static void fill_dlp_l_collocation_near_laplacebem3d(const uint *ridx,
													 const uint *cidx, pcbem3d bem,
													 bool ntrans, pamatrix N)
{
	pcsurface3d gr = bem->gr;
	real(*gr_x)[3] = (real(*)[3])gr->x;
	uint(*gr_t)[3] = (uint(*)[3])gr->t;
	plistnode *v2t = bem->v2t;
	field *aa = N->a;
	uint rows = ntrans ? N->cols : N->rows;
	uint cols = ntrans ? N->rows : N->cols;
	longindex ld = N->ld;
#if !defined(NDEBUG)
	uint vertices = gr->vertices;
	uint triangles = gr->triangles;
#endif
	ptri_list tl_c;

#ifdef USE_OPENMP
#pragma omp parallel if (!omp_in_parallel() && (cols > 64 * (1 << max_pardepth))) num_threads(1 << max_pardepth)
	{
#endif

		ptri_list tl1_c;
		pvert_list vl_c;
		plistnode v;
		uint j, jj, i, ii, l, cj, vv;
		uint *tri_j;
		real res[3];

#ifdef USE_OPENMP
#pragma omp single
		{
#endif
			clear_amatrix(N);
			tl_c = NULL;
			cj = 0;
			for (j = 0; j < cols; ++j)
			{
				jj = (cidx == NULL ? j : cidx[j]);
				for (v = v2t[jj], vv = v->data; v->next != NULL;
					 v = v->next, vv = v->data)
				{

					tl1_c = tl_c;
					while (tl1_c && tl1_c->t != vv)
					{
						tl1_c = tl1_c->next;
					}

					if (tl1_c == NULL)
					{
						tl1_c = tl_c = new_tri_list(tl_c);
						tl_c->t = vv;
						cj++;
					}

					tl1_c->vl = new_vert_list(tl1_c->vl);
					tl1_c->vl->v = j;
				}
			}

#ifdef USE_OPENMP
		}
#pragma omp for
#endif
		for (i = 0; i < rows; ++i)
		{
			//    	printf("completed: %d %\r", (int)((double)i / (double)rows * 100)); // not suitable for H-matrix calculation
			ii = (ridx == NULL ? i : ridx[i]);
			assert(ii < vertices);
			for (tl1_c = tl_c; tl1_c != NULL; tl1_c = tl1_c->next)
			{
				jj = tl1_c->t;
				assert(jj < triangles);

				tri_j = gr_t[jj];

				LindholmCWrapper(gr_x[ii], gr_x[tri_j[0]], gr_x[tri_j[1]],
								 gr_x[tri_j[2]], res);

				vl_c = tl1_c->vl;
				while (vl_c)
				{
					j = vl_c->v;
					if (j < cols)
					{
						jj = (cidx == NULL ? j : cidx[j]);
						for (l = 0; l < 3; ++l)
						{
							if (jj == tri_j[l])
							{
								if (ntrans)
								{
									aa[j + i * ld] += res[l];
								}
								else
								{
									aa[i + j * ld] += res[l];
								}
							}
						}
					}
					vl_c = vl_c->next;
				}
			}
		}

#ifdef USE_OPENMP
	}
#endif
	del_tri_list(tl_c);
}

void fill_nearfield_collocation_bem3d(pcbem3d bem, pamatrix A)
{
	bem->kernels->kernel_col(NULL, bem->gr->x, bem, true, A);
}

void assemble_fundamental_collocation_row_bem3d(const uint *idx,
												const real (*Z)[3], pcbem3d bem, bool trans, pamatrix A)
{
	real(*X)[3] = (real(*)[3])allocreal(3 * A->rows);
	uint i, ii;
	(void)trans;
	for (i = 0; i < A->rows; ++i)
	{
		ii = (idx != NULL ? idx[i] : i);
		X[i][0] = bem->gr->x[ii][0];
		X[i][1] = bem->gr->x[ii][1];
		X[i][2] = bem->gr->x[ii][2];
	}
	bem->kernels->fundamental(bem, X, Z, A);
	freemem(X);
}

void assemble_dnz_fundamental_collocation_row_bem3d(const uint *idx,
													const real (*Z)[3], const real (*N)[3], pcbem3d bem, bool trans, pamatrix A)
{
	(void)trans;
	real(*X)[3] = (real(*)[3])allocreal(3 * A->rows);
	uint i, ii;

	for (i = 0; i < A->rows; ++i)
	{
		ii = (idx != NULL ? idx[i] : i);
		X[i][0] = bem->gr->x[ii][0];
		X[i][1] = bem->gr->x[ii][1];
		X[i][2] = bem->gr->x[ii][2];
	}
	bem->kernels->dny_fundamental(bem, X, Z, N, A);
	freemem(X);
}

void assemble_lagrange_collocation_row_bem3d(const uint *idx, pcrealavector px,
											 pcrealavector py, pcrealavector pz, pcbem3d bem, pamatrix A)
{

	real(*X)[3] = (real(*)[3])allocreal(3 * A->rows);
	uint i, ii;

	for (i = 0; i < A->rows; ++i)
	{
		ii = (idx != NULL ? idx[i] : i);
		X[i][0] = bem->gr->x[ii][0];
		X[i][1] = bem->gr->x[ii][1];
		X[i][2] = bem->gr->x[ii][2];
	}

	assemble_bem3d_lagrange_amatrix(X, px, py, pz, bem, A);

	freemem(X);
}

pbem3d new_dlp_collocation_laplace_bem3d(pcsurface3d gr, uint q_regular,
										 uint q_singular,
										 basisfunctionbem3d basis_neumann,
										 basisfunctionbem3d basis_dirichlet)
{
	pbem3d bem;
	bem = new_dlp_laplace_bem3d(gr, q_regular, q_singular, basis_neumann,
								basis_dirichlet, 0.5);
	assert(basis_neumann == BASIS_LINEAR_BEM3D && basis_dirichlet == BASIS_LINEAR_BEM3D);
	bem->nearfield = fill_dlp_l_collocation_near_laplacebem3d;
	bem->nearfield_far = fill_dlp_l_collocation_near_laplacebem3d;
	bem->kernels->fundamental_row = assemble_fundamental_collocation_row_bem3d;
	bem->kernels->dnz_fundamental_row = assemble_dnz_fundamental_collocation_row_bem3d;
	bem->kernels->lagrange_row = assemble_lagrange_collocation_row_bem3d;
	return bem;
}

ph2matrix setup_H2coll_from_grid(pcsurface3d gr, uint q_reg, uint clf, real eta)
{
	pstopwatch sw;
	pbem3d bem;
	uint q_sing;
	pcluster root;
	pblock broot;
	phmatrix Kh;
	ph2matrix Kh2;
	real eps_aca;
	real eps_recomp;
	uint vertices;
	ptruncmode tm;
	tm = new_releucl_truncmode();
	q_sing = q_reg + 2;	
	eps_aca = 1.0e-8;	
	eps_recomp = 1.0e-5;

	sw = new_stopwatch();
	//	printf("Geometry has %d vertices, %d edges, %d triangles\n", gr->vertices,	gr->edges, gr->triangles);
	//	printf("================================\n");
	vertices = gr->vertices;

	bem = new_dlp_collocation_laplace_bem3d(gr, q_reg, q_sing, BASIS_LINEAR_BEM3D, BASIS_LINEAR_BEM3D);

	printf("Size of uncompressed matrix: %.3f MB.\n", sizeof(double) * vertices * vertices / 1024. / 1024.);
	printf("Calculating compressed matrix, part 1 of 2... ");
	fflush(stdout);
	start_stopwatch(sw);
	root = build_bem3d_cluster(bem, clf, BASIS_LINEAR_BEM3D);
	broot = build_strict_block(root, root, &eta, admissible_2_cluster);
	Kh = build_from_block_hmatrix(broot, 0);

	setup_hmatrix_aprx_hca_bem3d(bem, root, root, broot, 5, eps_aca);

	assemble_bem3d_hmatrix(bem, broot, Kh);
	printf("done.  (%.2f s)\n", stop_stopwatch(sw));
	printf("Size of H matrix: %.3f MB\n", getsize_hmatrix(Kh) / 1024.0 / 1024.0);

	//	printf("Convert H-matrix Kh to H2-matrix kh2:\n");
	printf("Continuing compression, part 2 of 2... ");
	fflush(stdout);
	start_stopwatch(sw);

	Kh2 = compress_hmatrix_h2matrix(Kh, tm, eps_recomp);
	printf("done.  (%.2f s)\n", stop_stopwatch(sw));
	printf("Size of H2 matrix: %.3f MB\n", (getsize_h2matrix(Kh2) + getsize_clusterbasis(Kh2->rb) + getsize_clusterbasis(Kh2->cb)) / 1024.0 / 1024.0);

	del_laplace_bem3d(bem);
	del_hmatrix(Kh);
	del_stopwatch(sw);

	printf("Writing H2 matrix to file %s\n", filenameC);
	write_cdfcomplete_h2matrix(Kh2, filenameC);
	return Kh2;
}
// ###################################### end collocation ###########################################

void setupH2MatrixC(double *diagData, int dim)
{
	uint q_reg, clf;
	real eta;
	uint bnx;
	char **a;
	init_h2lib(0, &a);
	q_reg = 3;
	eta = 2.0;
	clf = 32;

	if (readingH2)
	{
		printf("reading H2 matrix from file %s\n", filenameC);
		Kh2 = read_cdfcomplete_h2matrix(filenameC);
	}
	else
	{
		Kh2 = setup_H2coll_from_grid(gr2, q_reg, clf, eta);
	}
	bnx = getrows_h2matrix(Kh2);
	assert(bnx == getcols_h2matrix(Kh2));
	assert(dim == (int)bnx);
	diagPart = new_avector(dim);
	memcpy(diagPart->v, diagData, dim * sizeof(double));
	b = new_avector(bnx);
}

void deleteHmatrices()
{
	uint act_mat, act_vec;
	del_avector(b);
	del_h2matrix(Kh2);
	del_avector(diagPart);
	del_surface3d(gr2);
	act_mat = getactives_amatrix();
	act_vec = getactives_avector();
	if (act_mat + act_vec)
	{
		printf("-------- Memory Leak Detected --------------\n"
			   "  %u matrices and\n"
			   "  %u vectors still active\n",
			   act_mat, act_vec);
		printf("--------------------------------------------\n");
	}
	uninit_h2lib();
}
