#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

using std::istream;
using std::ostream;
using std::vector;

namespace zich
{

class Matrix 
{
    public:

        /* constructor */
        Matrix() = default;
        Matrix(int row, int col); // for zero mat
        Matrix (const std::vector<double>& mat, int row, int col);
        ~Matrix();

        // getters
        size_t get_row() const{ return this->row; }
        size_t get_col() const{ return this->col; }
        std::vector<std::vector<double>> get_mat() const { return this->my_mat; }
        double get_sum() const { return this->sum; } // sum all fields
        
        // setter
        void set_sum(double val){ this->sum = val; } // change the total sum
        void set_field(double element, unsigned int row, unsigned int col){ this->my_mat[row][col] = element; }

        /* arithmetic operators */
        Matrix operator-();
        Matrix operator-(const Matrix& mat_sub) const;
        Matrix& operator-=(const Matrix& mat_sub);

        Matrix operator+();
        Matrix operator+(const Matrix& mat_add) const;
        Matrix& operator+=(const Matrix& mat_add);

        Matrix operator*(const Matrix& mat_mul) const;
        Matrix& operator*=(const Matrix&);
        friend  Matrix operator*(const double num, const Matrix& mat_mul);
        Matrix operator*(double num) const;
        Matrix& operator*=(double num);

        Matrix& operator++();
        Matrix operator++(int dummy_flag_for_postfix_increament);

        Matrix& operator--();
        Matrix operator--(int dummy_flag_for_postfix_increament);

        /* boolean operators */
        bool operator==(const Matrix& mat) const;
        bool operator!=(const Matrix& mat) const;

        bool operator<(const Matrix& mat1) const;
        bool operator<=(const Matrix& mat1) const;

        bool operator>(const Matrix& mat1) const;
        bool operator>=(const Matrix& mat1) const;

        // globals
        friend istream& operator>>(istream &in, Matrix &mat);
        friend ostream& operator<<(ostream &out, const Matrix &mat);

    private:
        
        size_t row, col;
        double sum; // of all fields togheter
        std::vector<std::vector<double>> my_mat;

};        

}   