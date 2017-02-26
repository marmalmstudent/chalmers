#ifndef HANDIN5_H
#define HANDIN5_H

#include <vector>
#include <math.h>
#include <complex>
#include <armadillo>
#include <chrono>

static double PI = 3.141592653589793238463;
static std::complex<double> j(0, 1);
class Handin5
{
    public:
        Handin5();
        ~Handin5();
        void initTask1();
        void initTask2();
        void initTask3();
        void initTask4();
        void initTask5();
        void initTask6();
        void borderTX(std::complex<double> nStart,
                      std::complex<double> nAntiRefl,
                      std::vector<std::complex<double>>* out);
        void propTX(std::complex<double> layerThick,
                    std::complex<double> nLayer,
                    std::vector<std::complex<double>>* out);
        void dot33(const int nRows, const int nCols, const int nSlices,
                   std::vector<std::complex<double>>* mat1,
                   std::vector<std::complex<double>>* mat2,
                   std::vector<std::complex<double>>* out);
        void print(const int nRows, const int nCols, const int nSlices);


    protected:

    private:
        int _taskNbr;
        std::complex<double> _lbdaZero, _lbdaMax, _lbdaMin;
        std::vector<std::complex<double>> _lbda_0, _matrix;
};

#endif // HANDIN5_H
