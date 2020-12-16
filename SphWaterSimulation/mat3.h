#ifndef MAT3_H_
#define MAT3_H_

#include <vector>
#include <cassert>

#include "vec3.h"

class Mat3{
private:
    std::vector<Vec3> m_matrix;

public:
    Mat3(){
        m_matrix = std::vector<Vec3>{ Vec3(1.0, 0.0, 0.0), Vec3(0.0, 1.0, 0.0), Vec3(0.0, 0.0, 1.0) };
    }
    Mat3(const Mat3& a) = default;
    Mat3(std::vector<Vec3> other_mat){
        m_matrix = other_mat;
    }

    Mat3 Identity() const{
        return Mat3(std::vector<Vec3>{Vec3(1.0, 0.0, 0.0), Vec3(0.0, 1.0, 0.0), Vec3(0.0, 0.0, 1.0)});
    }
    Mat3 Zero() const{
        return Mat3(std::vector<Vec3>(3, Vec3(0.0, 0.0, 0.0)));
    }

    Vec3 multiply(Vec3 other_vector) const{
        return Vec3(m_matrix[0].dot(other_vector),
            m_matrix[1].dot(other_vector),
            m_matrix[2].dot(other_vector));
    }
    Mat3 T() const{
        std::vector<Vec3> new_matrix(3, Vec3());
        for(int i = 0; i < 3; i++){
            new_matrix[i] = this->col(i);
        }
        return Mat3(new_matrix);
    }
    Mat3 multiply(Mat3 other_mat) const{
        std::vector<Vec3> new_matrix(3, Vec3());
        for(int i = 0; i < 3; i++){
            new_matrix[i] = Vec3(m_matrix[i].dot(other_mat.col(0)),
                m_matrix[i].dot(other_mat.col(1)),
                m_matrix[i].dot(other_mat.col(2)));
        }
        return Mat3(new_matrix);
    }

    Mat3 operator*(const double a) const{
        std::vector<Vec3> new_matrix(3, Vec3());
        for(int i = 0; i < 3; i++){
            new_matrix[i] = m_matrix[i] * a;
        }
        return Mat3(new_matrix);
    }
    Mat3 operator+(const Mat3 a)const{
        std::vector<Vec3> new_matrix(3, Vec3());
        for(int i = 0; i < 3; i++){
            new_matrix[i] = m_matrix[i] + a.row(i);
        }
        return Mat3(new_matrix);
    }
    Mat3& operator+=(const Mat3& a){
        return *this = *this + a;
    }

    Vec3 row(int i) const{
        assert(i >= 0 && i < 3 && "Row index out of bounds!");
        return m_matrix[i];
    }

    Vec3 col(int j) const{
        assert(j >= 0 && j < 3 && "Col index out of bounds!");
        switch(j){
            case 0:
                return Vec3(m_matrix[0].getx(), m_matrix[1].getx(), m_matrix[2].getx());
            case 1:
                return Vec3(m_matrix[0].gety(), m_matrix[1].gety(), m_matrix[2].gety());
            case 2:
                return Vec3(m_matrix[0].getz(), m_matrix[1].getz(), m_matrix[2].getz());
        }
    }
};

#endif