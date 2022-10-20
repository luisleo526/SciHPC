/**
 * This file is part of the SparseMatrix library
 *
 * @license  MIT
 * @author   Petr Kessler (https://kesspess.cz)
 * @link     https://github.com/uestla/Sparse-Matrix
 */

#include <vector>
#include <iostream>
#include "exceptions.h"
#include "SparseMatrix.h"

// === CREATION ==============================================

SparseMatrix::SparseMatrix(int n) {
    this->construct(n, n);
}


SparseMatrix::SparseMatrix(int rows, int columns) {
    this->construct(rows, columns);
}


SparseMatrix::SparseMatrix(SparseMatrix &matrix) {
    this->deepCopy(matrix);
}


SparseMatrix &SparseMatrix::operator=(SparseMatrix &matrix) {
    if (&matrix != this) {
        this->destruct();
        this->deepCopy(matrix);
    }

    return *this;
}


void SparseMatrix::deepCopy(const SparseMatrix &matrix) {
    this->m = matrix.m;
    this->n = matrix.n;
    this->rows = new std::vector<int>(*(matrix.rows));

    if (matrix.vals != NULL) {
        this->cols = new std::vector<int>(*(matrix.cols));
        this->vals = new std::vector(*(matrix.vals));
    }
}


SparseMatrix::~SparseMatrix(void) {
    this->destruct();
}


void SparseMatrix::construct(int rows, int columns) {
    if (rows < 1 || columns < 1) {
        throw InvalidDimensionsException("Matrix dimensions cannot be zero or negative.");
    }

    this->m = rows;
    this->n = columns;

    this->vals = NULL;
    this->cols = NULL;
    this->rows = new std::vector<int>(rows + 1, 1);
    this->buffer = new DataType [columns];
}


void SparseMatrix::destruct() {
    if (this->vals != NULL) {
        delete this->vals;
        delete this->cols;
    }

    delete this->rows;
    delete this->buffer;
}


// === GETTERS / SETTERS ==============================================

int SparseMatrix::getRowCount() const {
    return this->m;
}


int SparseMatrix::getColumnCount() const {
    return this->n;
}


// === VALUES ==============================================


DataType SparseMatrix::get(int row, int col) const {
    this->validateCoordinates(row, col);

    int currCol;

    for (int pos = (*(this->rows))[row - 1] - 1; pos < (*(this->rows))[row] - 1; ++pos) {
        currCol = (*(this->cols))[pos];

        if (currCol == col) {
            return (*(this->vals))[pos];

        } else if (currCol > col) {
            break;
        }
    }

    return DataType();
}


SparseMatrix &SparseMatrix::set(DataType val, int row, int col) {
    this->validateCoordinates(row, col);

    int pos = (*(this->rows))[row - 1] - 1;
    int currCol = 0;

    for (; pos < (*(this->rows))[row] - 1; pos++) {
        currCol = (*(this->cols))[pos];

        if (currCol >= col) {
            break;
        }
    }

    if (currCol != col) {
        if (val != DataType()) {
            this->insert(pos, row, col, val);
        }

    } else if (val == DataType ()) {
        this->remove(pos, row);

    } else {
        (*(this->vals))[pos] = val;
    }

    return *this;
}


// === OPERATIONS ==============================================

void SparseMatrix::multiply(const DataType *x) const {

#pragma omp parallel for default(shared)
    for (int i = 0; i < this->m; i++) {
        DataType sum = 0.0;
        for (int j = (*(this->rows))[i]; j < (*(this->rows))[i + 1]; j++) {
            sum = sum + (*(this->vals))[j - 1] * x[(*(this->cols))[j - 1] - 1];
        }
        this->buffer[i] = sum;
    }
}

void SparseMatrix::operator*(DataType *x) const {
    this->multiply(x);
}


// === HELPERS / VALIDATORS ==============================================

void SparseMatrix::validateCoordinates(int row, int col) const {
    if (row < 1 || col < 1 || row > this->m || col > this->n) {
        throw InvalidCoordinatesException("Coordinates out of range.");
    }
}


void SparseMatrix::insert(int index, int row, int col, DataType val) {
    if (this->vals == NULL) {
        this->vals = new std::vector<DataType>(1, val);
        this->cols = new std::vector<int>(1, col);

    } else {
        this->vals->insert(this->vals->begin() + index, val);
        this->cols->insert(this->cols->begin() + index, col);
    }

    for (int i = row; i <= this->m; i++) {
        (*(this->rows))[i] += 1;
    }
}


void SparseMatrix::remove(int index, int row) {
    this->vals->erase(this->vals->begin() + index);
    this->cols->erase(this->cols->begin() + index);

    for (int i = row; i <= this->m; i++) {
        (*(this->rows))[i] -= 1;
    }
}


// === FRIEND FUNCTIONS =========================================

bool operator==(const SparseMatrix &a, const SparseMatrix &b) {
    return ((a.vals == NULL && b.vals == NULL)
            || (a.vals != NULL && b.vals != NULL && *(a.vals) == *(b.vals)))
           && ((a.cols == NULL && b.cols == NULL)
               || (a.cols != NULL && b.cols != NULL && *(a.cols) == *(b.cols)))
           && *(a.rows) == *(b.rows);
}


bool operator!=(const SparseMatrix &a, const SparseMatrix &b) {
    return !(a == b);
}


std::ostream &operator<<(std::ostream &os, const SparseMatrix &matrix) {
    for (int i = 1; i <= matrix.m; i++) {
        for (int j = 1; j <= matrix.n; j++) {
            if (j != 1) {
                os << " ";
            }

            os << matrix.get(i, j);
        }

        if (i < matrix.m) {
            os << std::endl;
        }
    }

    return os;
}