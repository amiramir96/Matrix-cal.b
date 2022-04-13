#include "Matrix.hpp"
#include <math.h>

constexpr double EPS = 1e-16; // will be used for double 

using std::istream;
using std::ostream;

namespace zich
{

    bool compare_proportion(const Matrix& mat1, const Matrix& mat2)
    /* this function check if col, row of two matrixes is completely the same! */
    { return (mat1.get_row() == mat2.get_row() && mat1.get_col() == mat2.get_col()); }

    bool possible_to_mul(const Matrix& mat1, const Matrix& mat2)
    /* this function verify that we can do the multiplication operation of 'mat1*mat2' via mat1_col == mat2_row */
    { return mat1.get_col() == mat2.get_row(); }

    Matrix::Matrix (const std::vector<double>& mat, int row, int col) 
    {        /* constructor, will make a 2dim vector filled by the mat input values */

        if (row <= 0 || col <= 0)
        {
            throw std::invalid_argument("cannot get negative or zero as sizes input\n");
        }
        if (row*col != mat.size())
        {
            throw std::invalid_argument("size of vector not equal to 2dim matrix\n");
        }

        double fields_value = 0;
        size_t temp = 0;
        for (size_t i=0; i < row; i++)
        {
            this->my_mat.push_back(std::vector<double>((size_t)col));
            for (size_t j=0; j < col; j++)
            {
                temp = (size_t)(i*(size_t)col +j);
                this->my_mat[i][j] = (mat[temp]); // the row idx i -> i rows, ea row 'col' fields, + j fields to curr field
                fields_value += mat[temp];
            }            
        }
        this->col = (size_t)col;
        this->row = (size_t)row;
        this->sum = fields_value;
    
    }

    Matrix::Matrix (int row, int col)
    { /* returns zero mat with row*col sizes */
        if (row <= 0 || col <= 0)
        {
            throw std::invalid_argument("cannot get negative or zero as sizes input\n");
        }

        for (int i=0; i < row; i++)
        { // create the 2nd dim
            this->my_mat.push_back(std::vector<double>((size_t)col));      
        }
        this->col = (size_t)col;
        this->row = (size_t)row;
        this->sum = 0;
    }

    Matrix::~Matrix()
    { /* destructor, will clear the 2dim vector */
        this->my_mat.clear();
    }

    /* arithmetic operators */
    Matrix Matrix::operator-(const Matrix& mat_sub) const
    { /* sub between two matrixes, set ans at new matrix same sizes (named res) */
        if (!compare_proportion(*this, mat_sub))
        {
            throw std::invalid_argument("matrixes sizes is not the same, cant support this operation\n");
        }

        Matrix res(get_row(), get_col()); // result mat same sizes
        for (size_t i=0; i < res.get_row(); i++) // minus cal
        {
            for (size_t j=0; j < res.get_col(); j++)
            {
                res.set_field(get_mat()[i][j] - mat_sub.get_mat()[i][j], i, j); // set ans at res matrix
            }
        }
        res.set_sum(get_sum() - mat_sub.get_sum()); // update totalsum field
        return res;
    }

    Matrix& Matrix::operator-=(const Matrix& mat_sub)
    { /* sub between two matrixes, origin matrix is the one being subbed */
        if (!compare_proportion(*this, mat_sub))
        {
            throw std::invalid_argument("matrixes sizes is not the same, cant support this operation\n");
        }

        for (size_t i=0; i < get_row(); i++) // minus cal
        {
            for (size_t j=0; j < get_col(); j++)
            {
                set_field(get_mat()[i][j] - mat_sub.get_mat()[i][j], i, j); // set ans at origin matrix
            }
        }
        set_sum(get_sum() - mat_sub.get_sum()); // update totalsum field
        return *this;
    }

    Matrix Matrix::operator-(){ /* return the negative matrix */ return (*this) * -1; }

    Matrix Matrix::operator+(){ /* just use mul func to return positive same mat */ return (*this)*1; }
    
    Matrix Matrix::operator+(const Matrix& mat_add) const
    { /* add between two matrixes, set ans at new matrix same sizes (named res) */
        if (!compare_proportion(*this, mat_add))
        {
            throw std::invalid_argument("matrixes sizes is not the same, cant support this operation\n");
        }

        Matrix res(get_row(), get_col()); // result mat same sizes
        for (size_t i=0; i < res.get_row(); i++) // additive cal
        {
            for (size_t j=0; j < res.get_col(); j++)
            {
                res.set_field(get_mat()[i][j] + mat_add.get_mat()[i][j], i, j); // set ans at res matrix
            }
        }
        res.set_sum(get_sum() + mat_add.get_sum()); // update totalsum field
        return res;
    }

    Matrix& Matrix::operator+=(const Matrix& mat_add)
    { /* add between two matrixes, origin matrix is the one being added */
        if (!compare_proportion(*this, mat_add))
        {
            throw std::invalid_argument("matrixes sizes is not the same, cant support this operation\n");
        }

        for (size_t i=0; i < get_row(); i++) // additive cal
        {
            for (size_t j=0; j < get_col(); j++)
            {
                set_field(get_mat()[i][j] + mat_add.get_mat()[i][j], i, j); // set ans at origin matrix
            }
        }
        set_sum(get_sum() + mat_add.get_sum()); // update totalsum field
        return *this;
    }
        
    Matrix Matrix::operator*(const Matrix& mat_mul) const
    { /* read here for definition of matrix mul: https://en.wikipedia.org/wiki/Matrix_multiplication
        the funct first verify that we can operate the mul 
        then, exec multiplication and save the values of ea field at result Matrix */
        if (!possible_to_mul(*this, mat_mul))
        {
            throw std::invalid_argument("mat1.col != mat2.row => cannot mul them! \n");
        }

        double temp_sum = 0;
        Matrix res(get_row(), mat_mul.get_col());
        for (size_t i=0; i < get_row(); i++)
        { // row of mat1
            for (size_t j=0; j < mat_mul.get_col(); j++)
            { // col of mat2  
                temp_sum = 0;
                for (size_t k=0; k < get_col(); k++)
                { // shared size (col mat1 // row mat 2)
                // current field of rest mat == [row of mat1][col of mat2]
                    temp_sum += get_mat()[i][k] * mat_mul.get_mat()[k][j];
                }
                res.set_field(temp_sum, i, j);
                res.set_sum(res.get_sum() + temp_sum); // add curr item val to the total sum of res (reminder - res created with sum = 0)
            }
        }
        return res;
    }
    Matrix& Matrix::operator*=(const Matrix& mat_mul)
    { /* will use mul bet 2 matrixes and then switch the origin (a.k.a '*this' ) to be res matrix */    
        Matrix res = (*this) * mat_mul;
        *this = res;
        return *this;
    }

    Matrix Matrix::operator*(double num) const
    { /* mul by scalar-'num' all fields of mat and set it at new matrix and returns it */
        Matrix res(get_row(), get_col()); // result matrix
        for (size_t i=0; i < res.get_row(); i++)
        {
            for (size_t j=0; j < res.get_col(); j++)
            {
                res.set_field(get_mat()[i][j] * num, i, j); // mul by scalar-'num'
            }
        }
        res.set_sum(get_sum()*num); // update total sum field
        return res;
    }

    Matrix operator*(double num, const Matrix& mat) 
    /* same as the above function, multiplie scalar and matrix is commutative operating ! */
    {
        Matrix res(mat.get_row(), mat.get_col()); // result matrix
        for (size_t i=0; i < res.get_row(); i++)
        {
            for (size_t j=0; j < res.get_col(); j++)
            {
                res.set_field(mat.get_mat()[i][j] * num, i, j); // mul by scalar-'num'
            }
        }
        res.set_sum(mat.get_sum()*num); // update total sum field
        return res;
    }

    Matrix& Matrix::operator*=(double num)
    { /* mul by scalar-'num' all fields of mat */
        for (size_t i=0; i < get_row(); i++)
        {
            for (size_t j=0; j < get_col(); j++)
            {
                set_field(get_mat()[i][j] * num, i, j); // mul by scalar-'num'
            }
        }
        set_sum(get_sum()*num); // update total sum field
        return *this;
    }

    Matrix& Matrix::operator++()
    { /* looks like ++mat, increase all fields of mat by 1 */
        std::vector<double> plus_one;
        for (size_t i=0; i < get_row(); i++)
        {
            for (size_t j=0; j < get_col(); j++)
            {
                set_field(get_mat()[i][j]+1, i, j);
            }
        }
        set_sum(get_sum() + (double) (get_row() * get_col())); // update sum
        return *this;
    }

    Matrix Matrix::operator++(int dummy_flag_for_postfix_increament)
    { /* looks like:  mat++
        http://www.cs.technion.ac.il/users/yechiel/c++-faq/increment-pre-post-overloading.html 
        implemented the function via this guide */
        Matrix temp = (*this);
        ++(*this);
        return temp;
    }

    Matrix& Matrix::operator--()
    { /* looks like --mat, decrease all fields of mat by 1 */
        for (size_t i=0; i < get_row(); i++)
        {
            for (size_t j=0; j < get_col(); j++)
            {
                set_field(get_mat()[i][j]-1, i, j);
            }
        }
        set_sum(get_sum() - (double) (get_row() * get_col())); // update sum
        return *this;
    }

    Matrix Matrix::operator--(int dummy_flag_for_postfix_increament)
    { /* looks like: mat--
        http://www.cs.technion.ac.il/users/yechiel/c++-faq/increment-pre-post-overloading.html 
        implemented the function via this guide */
        Matrix temp = *this;
        --(*this);
        return temp;
    }

    /* boolean operators */
    bool Matrix::operator==(const Matrix& mat) const
    { /* compare all fields of mat1, mat2, if and only if completly the same -> returns true*/
        bool not_same = false;
        if (!compare_proportion(*this, mat))
        { // sizes not same
            throw std::invalid_argument("matrixes sizes is not the same, cant support this operation\n");
        }
        for (size_t i=0; i < get_row(); i++)
        {
            for (size_t j=0; j < get_col(); j++)
            {
                if (std::fabs(get_mat()[i][j] - mat.get_mat()[i][j]) > EPS) // since we work with double
                {
                    return not_same; // one field not same is enough to return false 
                }
            }
        }
        return !not_same; // no false = true
    }

    bool Matrix::operator!=(const Matrix& mat) const
    /* not == is same as !=, at least one field of matrixes is not the same will return true */
    { return !((*this) == mat); }

    bool Matrix::operator<(const Matrix& mat) const
    /* mat1 < mat2 means mat1.sum < mat2.sum AND also they have both same row and col */
    { 
        if (!compare_proportion(*this, mat))
        { // sizes not same
            throw std::invalid_argument("matrixes sizes is not the same, cant support this operation\n");
        }
        return get_sum() < mat.get_sum();
    }
    
    bool Matrix::operator<=(const Matrix& mat) const
    /* mat1 <= mat2 means: mat1.sum <= mat2.sum equivalent to  not(mat1.sum > mat2.sum) */
    { return !((*this) > mat); }

    bool Matrix::operator>(const Matrix& mat) const
    /* mat1 > mat2 means mat1.sum > mat2.sum AND also they have both same row and col */
    { 
        if (!compare_proportion(*this, mat))
        { // sizes not same
            throw std::invalid_argument("matrixes sizes is not the same, cant support this operation\n");
        }
        return get_sum() > mat.get_sum();
    }

    bool Matrix::operator>=(const Matrix& mat) const
    /* mat1 >= mat2 means: mat1.sum >= mat2.sum equivalent to  not(mat1.sum < mat2.sum) */
    { return !((*this) < mat); }

    // Global Operators
    istream& operator>>(istream &in, Matrix &mat) 
    { /* takes cin and construct it to Matrix to &mat 
        format of cin shall be: '[1 0 0], [0 1 0], [0 0 1]' fow ex<- identity matrix 3x3
        the responsibility is completly on the USER!
            if row*col != amount of elements in string the func will throw 'invalid_argument' exception! 
      */

      /*
        how its gonna work:
            1- hold the whole cin at string var
            2- split string via ', ' to rows
            3- remove from ea row the '[' <- strart, end -> ']'
            4- xfer ea string of row from format "num1 num2 num3 .... num_col" to double items and add it to vector<double>
            5- create the matrix
       */

        // set string at input_mat var
        std::string input_mat;
        std::getline(in, input_mat);
        // credit to 841 comment, helped alot! https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c
        // its a way to split string properly

        // --- block of spliting the whole string to rows! via ',' deli ---
        std::string delimeter = ",";
        size_t pos = 0;
        std::vector<std::string> tokens; // ea string represent a row of matrix formated "[num1 num2 .... num_col]"

        while((pos = input_mat.find(delimeter)) != std::string::npos)
        {
            tokens.push_back(input_mat.substr(0, pos));
            input_mat.erase(0, pos+2);
        }
        tokens.push_back(input_mat.substr(0, input_mat.size()));
        // ------ end of block ------

        int row = tokens.size(); // splited to rows, size of string vector is amount of rows

        // --- removes '[' ']' from start and end of ea string that represent a row
        std::vector<double> mat_vector;
        for (size_t j=0; j < row; j++)
        {
            tokens[j].erase(0, 1);
            tokens[j].erase(tokens[j].size()-1, 1);
            tokens[j].append(" ");
        }
        // ----- end of block -----

        // this block split again, this time each row to double elements via ' ' deli
        for (size_t i=0; i < row; i++)
        {
            while ((pos = tokens[i].find(' ')) != std::string::npos) 
            {
                mat_vector.push_back(std::stod(tokens[i].substr(0, pos)));
                tokens[i].erase(0, pos+1);
            }
        }
        // ----- end of block -----

        if (row == 0){throw std::invalid_argument("row is zero, invalid istream\n"); } // for the 'tidy' tester, without this line there is 'warning' that row can be zero :(

        int col = ((int)(mat_vector.size())) / row; // amount_of_cols = all_elements / amount_of_rows 
        mat = Matrix(mat_vector, row, col);
        return in; 
    }

    ostream& operator<<(ostream &out, const Matrix &mat)
    { /* creates string from a matrix for example: (identity matrix 3x3)
         [1 0 0]
         [0 1 0]
         [0 0 1]
        */

        std::string mat_str;
        std::string temp;

        // this loops shall create a string formated to matrix 
        for (size_t i=0; i < mat.get_row(); i++)
        {
            mat_str += "[";
            size_t j=0;
            for (;j < mat.get_col(); j++)
            {

                // this block will handle situations of 'negative zero' in cpp (so printing shall be more sympathic for user)
                if (mat.my_mat[i][j] == 0)
                {
                    temp = "0.0";
                }
                else 
                {
                    temp = std::to_string(mat.my_mat[i][j]);
                }
                // ----- end of block ---------

                // this block handle trail of zero in the end of the floating number
                while (temp.at(temp.size()-1) == '0')
                {
                    temp.erase(temp.size()-1);
                }
                if (temp.at(temp.size()-1) == '.')
                {
                    temp.erase(temp.size()-1);
                }
                //---------- end of block --------------------------------
                
                mat_str += temp;
                mat_str += " ";
            }
            if (j > 0)
            { // erase last whitespace if exist
                mat_str.erase(mat_str.size()-1);
            }
            mat_str += "]\n";
        }

        mat_str.erase(mat_str.size()-1); // remove last downline, means cursor will be on the end of the matrix printing and not line below!
        out << mat_str;
        return out;
    }
}