#include "dns.h"
#include "mpi.h"
#include <cmath>
#include <iostream>

using namespace std;

#define iDIM 0
#define jDIM 1
#define kDIM 2


void create_topology(MeshInfo *meshInfo) {
    int dims[3];
    int periods[3];
    dims[iDIM] = meshInfo->dim_size;
    dims[jDIM] = meshInfo->dim_size;
    dims[kDIM] = meshInfo->dim_size;
    periods[iDIM] = periods[jDIM] = periods[kDIM] = 1;
    meshInfo->mesh3d = MPI_COMM_WORLD;

    int err = MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0,
                              &meshInfo->mesh3d);
    if (err != 0) {
        cout << "Error in MPI_Cart_create" << endl;
        exit(1);
    }

    MPI_Cart_coords(meshInfo->mesh3d, meshInfo->my_rank, 3, meshInfo->coords);

    dims[iDIM] = dims[kDIM] = 1;
    dims[jDIM] = 0;

    MPI_Cart_sub(meshInfo->mesh3d, dims, &meshInfo->mesh_ik);

    dims[iDIM] = dims[kDIM] = 0;
    dims[jDIM] = 1;

    MPI_Cart_sub(meshInfo->mesh3d, dims, &meshInfo->ring_j);

    dims[jDIM] = dims[kDIM] = 1;
    dims[iDIM] = 0;

    MPI_Cart_sub(meshInfo->mesh3d, dims, &meshInfo->mesh_jk);

    dims[jDIM] = dims[kDIM] = 0;
    dims[iDIM] = 1;

    MPI_Cart_sub(meshInfo->mesh3d, dims, &meshInfo->ring_i);

    dims[iDIM] = dims[jDIM] = 1;
    dims[kDIM] = 0;

    MPI_Cart_sub(meshInfo->mesh3d, dims, &meshInfo->mesh_ij);

    dims[iDIM] = dims[jDIM] = 0;
    dims[kDIM] = 1;

    MPI_Cart_sub(meshInfo->mesh3d, dims, &meshInfo->ring_k);
}

Matrix *MultiplyMatrixesDNS(Matrix *A, Matrix *B, MeshInfo *t) {
    int *buf = NULL, *matrix = NULL;;
    int block_size;

    Matrix *Ablock = DistributeLeftMatrix(A, t);   //分发矩阵A, 每个节点分到的小矩阵放在Ablock中
    Matrix *Bblock = DistributeRightMatrix(B, t);   //分发矩阵B， 每个节点分到的矩阵放在Bblock中
    Matrix *Cblock = MultiplyMatrixes(Ablock, Bblock);   //Cblock = Ablock * Bblock
    Matrix *C = NULL;    //用于存放C矩阵

    block_size = t->matrix_size / t->dim_size;   //获取子块大小

    buf = (int *) malloc(sizeof(int) * Cblock->n * Cblock->n);   //分配缓存用于保存子块
    memset(buf, 0, sizeof(int) * Cblock->n * Cblock->n);

    //在J维度将数据规约
    MPI_Reduce(MatrixToArrByCols(Cblock), buf, Cblock->n * Cblock->n, MPI_INT, MPI_SUM, 0, t->ring_j);

    if (t->coords[jDIM] == 0) {
        matrix = (int *) malloc(sizeof(int) * t->matrix_size * t->matrix_size);

        //将I维度和K维度的数据保存在j维度为1的所有节点里
        MPI_Gather(buf, block_size * block_size, MPI_INT, matrix, block_size * block_size, MPI_INT, 0, t->mesh_ik);

        if (buf != NULL) free(buf);
    }

    if (t->my_rank == 0) {
        C = ArrBlocksToMatrix(matrix, block_size, block_size, t->matrix_size, t->matrix_size);  //将数组转成矩阵
        free(matrix);
    }


    FreeMatrix(Ablock);
    FreeMatrix(Bblock);
    FreeMatrix(Cblock);

    return C;
}

Matrix *DistributeLeftMatrix(Matrix *A, MeshInfo *topology) {
    Matrix *block;
    int blockSize = topology->matrix_size / topology->dim_size;       //块大小
    int *buf = (int *) malloc(sizeof(int) * blockSize * blockSize);   //为块分配空间
    int *send = NULL;

    if (topology->my_rank == 0) {
        send = MatrixToArrBlocksRows(A, blockSize, blockSize);        //将矩阵A转成线性数组
    }

    if (topology->coords[kDIM] == 0) {
        MPI_Scatter(send, blockSize * blockSize, MPI_INT, buf, blockSize * blockSize, MPI_INT, 0,
                    topology->mesh_ij);     //将矩阵A的数据分发到k维度为0的那一层
    }

    MPI_Bcast(buf, blockSize * blockSize, MPI_INT, 0, topology->ring_k);    //k维度为0的将数据广播到整个k维度
    block = ArrToMatrixByRows(buf, blockSize, blockSize);   //将数组转成矩阵

    //释放资源
    if (topology->my_rank == 0) {
        free(send);
    }

    free(buf);
    return block;

}

Matrix *DistributeRightMatrix(Matrix *A, MeshInfo *topology) {
    Matrix *block;
    int blockSize = topology->matrix_size / topology->dim_size;  //块大小
    int *buf = (int *) malloc(sizeof(int) * blockSize * blockSize);  //为块分配空间
    int *send = NULL;
    if (topology->my_rank == 0) {
        send = MatrixToArrBlocksRows(A, blockSize, blockSize);  //将矩阵A转成线性数组
    }

    if (topology->coords[iDIM] == 0) {
        MPI_Scatter(send, blockSize * blockSize, MPI_INT, buf, blockSize * blockSize, MPI_INT, 0,
                    topology->mesh_jk);    //将矩阵A的数据分发到k维度为0的那一层
    }


    MPI_Bcast(buf, blockSize * blockSize, MPI_INT, 0, topology->ring_i);  //k维度为0的将数据广播到整个k维度
    block = ArrToMatrixByRows(buf, blockSize, blockSize);  //将数组转成矩阵

    //释放资源
    if (topology->my_rank == 0) {
        free(send);
    }

    free(buf);
    return block;
}

int *MatrixToArrByRows(Matrix *matrix) {
    int *a = (int *) malloc(sizeof(int) * matrix->n * matrix->m);
    for (int i = 0; i < matrix->m; i++) {
        for (int j = 0; j < matrix->n; j++) {
            a[i * matrix->n + j] = matrix->data[i][j];
        }
    }
    return a;
}

int *MatrixToArrByCols(Matrix *matrix) {
    int *matrix_col;
    matrix_col = (int *) malloc(sizeof(int) * matrix->n * matrix->m);
    for (int i = 0; i < matrix->n; i++) {
        for (int j = 0; j < matrix->m; j++) {
            matrix_col[i * matrix->m + j] = matrix->data[j][i];
        }
    }

    return matrix_col;
}

Matrix *ArrToMatrixByRows(int *matrix_row, int n, int m) {
    Matrix *matrix = CreateMatrix(n, m);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            matrix->data[i][j] = matrix_row[i * n + j];
        }
    }
    return matrix;
}

Matrix *ArrToMatrixByCols(int *a, int n, int m) {
    Matrix *A = CreateMatrix(n, m);
    for (int i = 0; i < A->n; i++) {
        for (int j = 0; j < A->m; j++) {
            A->data[j][i] = a[i * A->m + j];
        }
    }
    return A;
}

int *MatrixToArrBlocksRows(Matrix *A, int blockrow, int blockcol) {
    int fullrow = A->m / blockrow;
    int fullcol = A->n / blockcol;
    int remainsrow = A->m % blockrow;
    int remainscol = A->n % blockcol;
    int *a = (int *) malloc(sizeof(int) * A->n * A->m);
    int count = 0;
    int i, j;
    for (i = 0; i < fullrow; i++) {
        for (j = 0; j < fullcol; j++) {
            CopyBlockToArr(A, a, i * blockrow, j * blockcol, blockrow, blockcol, count);
            count += blockrow * blockcol;
        }

        CopyBlockToArr(A, a, i * blockrow, blockcol * fullcol, blockrow, remainscol, count);
        count += blockrow * remainscol;


    }
    for (j = 0; j < fullcol; j++) {
        CopyBlockToArr(A, a, i * blockrow, j * blockcol, remainsrow, blockcol, count);
        count += remainsrow * blockcol;
    }
    CopyBlockToArr(A, a, fullcol * blockrow, blockcol * fullcol, remainsrow, remainscol, count);


    return a;
}

int *MatrixToArrBlocksCols(Matrix *A, int blockrow, int blockcol) {
    int fullrow = A->m / blockrow;
    int fullcol = A->n / blockcol;
    int remainsrow = A->m % blockrow;
    int remainscol = A->n % blockcol;
    int *a = (int *) malloc(sizeof(int) * A->n * A->m);
    int count = 0;
    int i, j;
    for (j = 0; j < fullcol; j++) {
        for (i = 0; i < fullrow; i++) {
            CopyBlockToArr(A, a, i * blockrow, j * blockcol, blockrow, blockcol, count);
            count += blockrow * blockcol;
        }

        CopyBlockToArr(A, a, i * blockrow, j * blockcol, remainsrow, blockcol, count);
        count += remainsrow * blockcol;
        CopyBlockToArr(A, a, i * blockrow, blockcol * fullcol, blockrow, remainscol, count);
        count += blockrow * remainscol;


    }
    for (i = 0; i < fullrow; i++) {
        CopyBlockToArr(A, a, i * blockrow, blockcol * fullcol, blockrow, remainscol, count);
        count += blockrow * remainscol;
    }
    CopyBlockToArr(A, a, fullcol * blockrow, blockcol * fullcol, remainsrow, remainscol, count);


    return a;
}

void CopyBlockToArr(Matrix *A, int *a, int row, int col, int height, int width, int position) {
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            a[position + i * width + j] = A->data[row + i][col + j];
        }
    }
}

void CopyArrBlocktoMatrix(Matrix *A, int *a, int row, int col, int height, int width, int position) {
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            A->data[col + j][row + i] = a[position + i * width + j];
        }
    }
}

Matrix *ArrBlocksToMatrix(int *a, int blockrow, int blockcol, int n, int m) {
    Matrix *A = CreateMatrix(n, m);
    int fullrow = A->m / blockrow;
    int fullcol = A->n / blockcol;
    int remainsrow = A->m % blockrow;
    int remainscol = A->n % blockcol;

    int count = 0;
    int i, j;
    for (j = 0; j < fullcol; j++) {
        for (i = 0; i < fullrow; i++) {
            CopyArrBlocktoMatrix(A, a, i * blockrow, j * blockcol, blockrow, blockcol, count);
            count += blockrow * blockcol;
        }

        CopyArrBlocktoMatrix(A, a, i * blockrow, j * blockcol, remainsrow, blockcol, count);
        count += remainsrow * blockcol;
        CopyArrBlocktoMatrix(A, a, i * blockrow, blockcol * fullcol, blockrow, remainscol, count);
        count += blockrow * remainscol;


    }
    for (i = 0; i < fullrow; i++) {
        CopyArrBlocktoMatrix(A, a, i * blockrow, blockcol * fullcol, blockrow, remainscol, count);
        count += blockrow * remainscol;
    }
    CopyArrBlocktoMatrix(A, a, fullcol * blockrow, blockcol * fullcol, remainsrow, remainscol, count);

    return A;
}



