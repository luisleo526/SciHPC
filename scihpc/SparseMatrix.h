/**
 * This file is part of the SparseMatrix library
 *
 * @license  MIT
 * @author   Petr Kessler (https://kesspess.cz)
 * @link     https://github.com/uestla/Sparse-Matrix
 */

#ifndef __SPARSEMATRIX_H__

#define    __SPARSEMATRIX_H__

#include <vector>
#include <iostream>
#include "global.h"

class SparseMatrix {

public:

    DataType *buffer;

    // === CREATION ==============================================

    explicit SparseMatrix(int n); // square matrix n×n
    SparseMatrix(int rows, int columns); // general matrix

    SparseMatrix(SparseMatrix &m); // copy constructor
    SparseMatrix &operator=(SparseMatrix &m);

    ~SparseMatrix();


    // === GETTERS / SETTERS ==============================================

    [[nodiscard]] int getRowCount() const;

    [[nodiscard]] int getColumnCount() const;


    // === VALUES ==============================================

    [[nodiscard]] DataType get(int row, int col) const;

    SparseMatrix &set(DataType val, int row, int col);


    // === OPERATIONS ==============================================

    void multiply(const DataType *x) const;

    void operator*(DataType *x) const;

    // === FRIEND FUNCTIONS =========================================

    friend bool operator==(const SparseMatrix &a, const SparseMatrix &b);

    friend bool operator!=(const SparseMatrix &a, const SparseMatrix &b);

    friend std::ostream &operator<<(std::ostream &os, const SparseMatrix &matrix);


protected:

    int m, n;

    std::vector<DataType> *vals;
    std::vector<int> *rows, *cols;


    // === HELPERS / VALIDATORS ==============================================

    void construct(int m, int n);

    void destruct();

    void deepCopy(const SparseMatrix &m);

    void validateCoordinates(int row, int col) const;

    void insert(int index, int row, int col, DataType val);

    void remove(int index, int row);

};


#endif
