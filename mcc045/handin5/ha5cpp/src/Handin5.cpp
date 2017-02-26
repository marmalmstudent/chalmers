#include "Handin5.h"

Handin5::Handin5()
{
    //ctor
}

Handin5::~Handin5()
{
    //dtor
}
void Handin5::borderTX(std::complex<double> nLeft,
                       std::complex<double> nRight,
                       std::vector<std::complex<double>>* out)
{
    out->push_back((nRight+nLeft) / (((std::complex<double>)2)*nRight));
    out->push_back((nRight-nLeft) / (((std::complex<double>)2)*nRight));
    out->push_back((nRight-nLeft) / (((std::complex<double>)2)*nRight));
    out->push_back((nRight+nLeft) / (((std::complex<double>)2)*nRight));
}
void Handin5::propTX(std::complex<double> layerThick,
                     std::complex<double> nLayer,
                     std::vector<std::complex<double>>* out)
{
    std::vector<std::complex<double>> kZero;
    std::vector<std::complex<double>>::iterator itr;
    for (itr = _lbda_0.begin(); itr != _lbda_0.end(); ++itr)
        kZero.push_back((std::complex<double>)(2*PI)/(*itr));

    // fill S(1,1) with e^(j*k_0*n*d)
    std::vector<std::complex<double>> m11;
    for (itr = kZero.begin(); itr != kZero.end(); ++itr)
        m11.push_back(std::exp(j*(*itr)*nLayer*layerThick));
    // fill (1,2) with 0
    std::vector<std::complex<double>> m12;
    for (itr = kZero.begin(); itr != kZero.end(); ++itr)
        m12.push_back((std::complex<double>)0);
    // fill (2,1) with 0
    std::vector<std::complex<double>> m21;
    for (itr = kZero.begin(); itr != kZero.end(); ++itr)
        m21.push_back(std::complex<double>(0, 0));
    // fill (2,2) with e^(-j*k_0*n*d)
    std::vector<std::complex<double>> m22;
    for (itr = kZero.begin(); itr != kZero.end(); ++itr)
        m22.push_back(std::exp(-j*(*itr)*nLayer*layerThick));

    out->insert(out->begin(), m11.begin(), m11.end());
    out->insert(out->end(), m12.begin(), m12.end());
    out->insert(out->end(), m21.begin(), m21.end());
    out->insert(out->end(), m22.begin(), m22.end());
}

void Handin5::dot33(const int nRows, const int nCols, const int nSlices,
                    std::vector<std::complex<double>>* mat1,
                    std::vector<std::complex<double>>* mat2,
                    std::vector<std::complex<double>>* out)
{
    int rowSkip = nCols*nSlices;
    int colSkip = nSlices;
    int currOutIdx = 0;
    for (int r = 0; r < nRows; ++r)
    {
        for (int c = 0; c < nCols; ++c)
        {
            for (int s = 0; s < nSlices; ++s)
            {
                currOutIdx = r*rowSkip + c*colSkip + s;
                out->push_back(0);
                for (int k = 0; k < nCols; ++k)
                {
                    //out->at(currOutIdx) += mat1->at(r*rowSkip + k*colSkip + s) * mat2->at(k*rowSkip + c*colSkip + s); // 33
                    //out->at(currOutIdx) += mat1->at(r*rowSkip + k*colSkip + s) * mat2->at(k*nCols + c); // 32
                    out->at(currOutIdx) += mat1->at(r*nCols + k) * mat2->at(k*rowSkip + c*colSkip + s); // 23
                    //out->at(currOutIdx) += mat1->at(r*nCols + k) * mat2->at(k*nCols + c); // 22
                }
            }
        }
    }

}
void Handin5::print(const int nRows, const int nCols, const int nSlices)
{
    int rowSkip = nCols*nSlices;
    int colSkip = nSlices;
    for (int r = 0; r < nRows; ++r)
    {
        for (int c = 0; c < nCols; ++c)
        {
            std::cout << "[";
            for (int s = 0; s < nSlices; ++s)
            {
                std::cout << std::real(_matrix[r*rowSkip + c*colSkip + s]) << " + j"
                    << std::imag(_matrix[r*rowSkip + c*colSkip + s]) << ", ";
            }
            std::cout << "]" << std::endl;
        }
        std::cout << std::endl;
    }
}
void Handin5::initTask1()
{
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    _taskNbr = 1;
    _lbdaZero = (std::complex<double>)550e-9;
    _lbdaMin = (std::complex<double>)400e-9;
    _lbdaMax = (std::complex<double>)700e-9;
    for (double i = std::real(_lbdaMin); i < std::real(_lbdaMax); i+=1e-10)
    {
        _lbda_0.push_back((std::complex<double>)i);
    }
    const int nSlices = _lbda_0.size();
    std::complex<double> nStart(1);
    std::complex<double> nSubst(1.5);
    std::complex<double> nAntiRefl = sqrt(nSubst);
    std::complex<double> dAntiRefl = _lbdaZero/(((std::complex<double>)4)*nAntiRefl);
    std::vector<std::complex<double>> mat1;
    std::vector<std::complex<double>> mat2;
    borderTX(nStart, nSubst, &mat1);
    propTX(dAntiRefl, nAntiRefl, &mat2);
    std::vector<std::complex<double>>::iterator itr;
    /*
    for (itr = mat1.begin(); itr != mat1.end(); ++itr)
    {
        std::cout << *itr << std::endl;
    }
    std::cout << std::endl;
    for (itr = mat2.begin(); itr != mat2.end(); ++itr)
    {
        std::cout << *itr << std::endl;
    }
    */
    dot33(2, 2, nSlices, &mat1, &mat2, &_matrix);
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds time_span = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1);
    std::cout << "C++ Time: " << time_span.count() << " ns" << std::endl;
    //print(2, 2, nSlices);
}
void Handin5::initTask2()
{

}
void Handin5::initTask3()
{

}
void Handin5::initTask4()
{

}
void Handin5::initTask5()
{

}
void Handin5::initTask6()
{

}
