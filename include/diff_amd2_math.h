#ifndef DIFF_AMD2_MATH_H
#define DIFF_AMD2_MATH_H

#include <cmath>
#include <string>

namespace AMD2_math
{
    inline std::string base32_token(const std::string& value)
    {
        static constexpr char alphabet[]="ABCDEFGHIJKLMNOPQRSTUVWXYZ234567";
        std::string out="s_"; unsigned int buffer=0,bits=0;
        for(const unsigned char ch:value){buffer=(buffer<<8)|ch;bits+=8;while(bits>=5){bits-=5;out.push_back(alphabet[(buffer>>bits)&31]);}}
        if(bits)out.push_back(alphabet[(buffer<<(5-bits))&31]);
        return out;
    }

    inline double shear_contraction(const double A[3][3],const double S[3][3],const double M[3])
    {
        double sum=0.;
        for (int k=0;k<3;++k)
            for (int i=0;i<3;++i)
                for (int j=0;j<3;++j)
                    sum+=M[k]*A[i][k]*A[j][k]*S[i][j];
        return sum;
    }

    inline double buoyancy_contraction(const double A[3][3],const double q[3],const double d[3],const double M[3])
    {
        double sum=0.;
        for (int k=0;k<3;++k)
        {
            double projection=0.;
            for (int i=0;i<3;++i) projection+=d[i]*A[i][k];
            sum+=M[k]*projection*q[k];
        }
        return sum;
    }

    inline double scalar_contraction(const double A[3][3],const double p[3],const double M[3])
    {
        double sum=0.;
        for (int k=0;k<3;++k)
            for (int i=0;i<3;++i)
                sum+=M[k]*A[i][k]*p[k]*p[i];
        return sum;
    }
}

#endif
