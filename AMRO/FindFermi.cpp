#include "stdafx.h"
#include "DataExtractor.h"
#include "FindFermi.h"
using namespace std;

//int FindFermi::func(Ipp64f * params, Ipp64f * argkz, Ipp64f * argCos, Ipp64f * argSin, Ipp64f * r, int length, Ipp64f * temp, Ipp64f * out)
int FindFermi::func(Ipp64f *params, Ipp64f * argkz, Ipp64f * kx, Ipp64f * ky, int length, Ipp64f *temp, Ipp64f *out) {

	ippsMulC_64f(kx, 3.747665940, temp, length);
	//ippsMul_64f_I(r, temp, length);
	vdCos(length, temp, &temp[1 * length]); // cos cos
	ippsMulC_64f(ky, 3.747665940, temp, length);
	//ippsMul_64f_I(r, temp, length);
	vdCos(length, temp, &temp[2 * length]); // cos sin
	ippsMulC_64f(kx, 3.747665940 / 2, temp, length);
	//ippsMul_64f_I(r, temp, length);
	vdCos(length, temp, &temp[3 * length]); // cos cos/2
	ippsMulC_64f(ky, 3.747665940 / 2, temp, length);
	//ippsMul_64f_I(r, temp, length);
	vdCos(length, temp, &temp[4 * length]); // cos sin/2
	ippsMulC_64f(kx, 3.747665940 * 2, temp, length);
	//ippsMul_64f_I(r, temp, length);
	vdCos(length, temp, &temp[5 * length]); // cos 2 cos
	ippsMulC_64f(ky, 3.747665940 * 2, temp, length);
	//ippsMul_64f_I(r, temp, length);
	vdCos(length, temp, &temp[6 * length]); // cos 2 sin
	ippsMulC_64f(argkz, 0.5, &temp[8 * length], length); //kzc / 2
	vdCos(length, &temp[8 * length], &temp[7 * length]); // cos kzc/2
	vdCos(length, &temp[8 * length], &temp[9 * length]); // cos kzc ***made same as above, cos kzc/2

	ippsAdd_64f(&temp[5 * length], &temp[6 * length], temp, length);// param 5
	ippsMulC_64f(temp, -35164.83516*params[5 - 1], out, length);

	ippsMul_64f(&temp[1 * length], &temp[2 * length], temp, length);// param 4
	ippsMulC_64f_I(-35164.83516 * 2 * params[4 - 1], temp, length);
	ippsAdd_64f_I(temp, out, length);

	ippsAdd_64f(&temp[1 * length], &temp[2 * length], temp, length);// param 3
	ippsMulC_64f_I(-35164.83516 * params[3 - 1], temp, length);
	ippsAdd_64f_I(temp, out, length);
	ippsAddC_64f_I(35164.83516 / 2 * params[2 - 1], out, length);// param 2
	ippsSub_64f(&temp[2 * length], &temp[1 * length], temp, length);// param 6
	ippsSqr_64f_I(temp, length); //square
	ippsMul_64f_I(&temp[3 * length], temp, length); // mult by cos cos/2
	ippsMul_64f_I(&temp[4 * length], temp, length); // mult by cos sin/2
	ippsMul_64f_I(&temp[7 * length], temp, length); // mult by cos  kz/2
	ippsMulC_64f_I(-35164.83516 * params[6 - 1], temp, length);
	ippsMulC_64f_I(-35164.83516 *params[7 - 1], &temp[9 * length], length);
	ippsAdd_64f_I(temp, out, length);
	ippsAdd_64f_I(&temp[9 * length], out, length);

	return 0;
}
/*
int FindFermi::interfunc(Ipp64f *theta, Ipp64f *kx, Ipp64f *ky, int Laylength, int clength, int flength, Ipp64f *temp, Ipp64f *temp2, Ipp64f *out) {
	int count;
	for (int i = 0; i < clength; ++i) {
		ippsSub_64f(&kx[i*Laylength], &kx[i*Laylength + 1], &temp[2 * 0 * (Laylength - 1)], Laylength - 1);//kx[i+1]-kx[i]
		ippsSub_64f(&ky[i*Laylength], &ky[i*Laylength + 1], &temp[(2 * 0 + 1) * (Laylength - 1)], Laylength - 1);//ky[i+1]-ky[i]
		ippsMul_64f_I(&temp[2 * 0 * (Laylength - 1)], &temp[2 * 0 * (Laylength)], Laylength - 1);//(kx[i+1]-kx[i])^2
		ippsMul_64f_I(&temp[(2 * 0 + 1) * (Laylength - 1)], &temp[(2 * 0 + 1) * (Laylength - 1)], Laylength - 1);//(ky[i+1]-ky[i])^2
		ippsAdd_64f_I(&temp[2 * 0 * (Laylength - 1)], &temp[(2 * 0 + 1) * (Laylength - 1)], Laylength - 1);//(kx[i+1]-kx[i])^2+(ky[i+1]-ky[i])^2
		ippsSqrt_64f_I(&temp[(2 * 0 + 1) * (Laylength - 1)], Laylength - 1);//Sqrt[(kx[i+1]-kx[i])^2+(ky[i+1]-ky[i])^2]
		for (int j = 0; j < Laylength; ++j) {
			temp[j] = 0;

			for (int k = 0; k < j; ++k) {
				temp[j] = temp[j] + temp[(2 * 0 + 1) * (Laylength - 1) + k];//cumulative sum
			}
		}
		for (int j = 0; j < flength; ++j) {
			temp2[j] = temp[Laylength - 1] / flength * j;
		}
		count = 1;
		out[i * flength] = 0;
		for (int k = 0; k < Laylength; ++k) {
			if (count > flength - 1) break;
			if (temp[k] >= temp2[count]) {
				out[i * flength + count] = theta[k - 1] + (theta[k] - theta[k - 1])*(temp2[count] - temp[k - 1]) / (temp[k] - temp[k - 1]);
				++count;

			}
		}




	}
}
*/
int FindFermi::interfunc( Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int *Laylength,int Laynum, Ipp64f distance, Ipp64f *temp, Ipp64f *temp2, Ipp64f *outx,Ipp64f *outy,Ipp64f*outz) {
	int count;
	int cumLaylength = 0;
	int finalPoint = 0;
	for (int i = 0; i < Laynum; ++i) {
		ippsSub_64f(&kx[cumLaylength], &kx[cumLaylength + 1], temp, Laylength[i] - 1);//kx[i+1]-kx[i]
		ippsSub_64f(&ky[cumLaylength], &ky[cumLaylength + 1], &temp[(Laylength[i] - 1)],  Laylength[i] - 1);//ky[i+1]-ky[i]
		ippsMul_64f_I(&temp[2 * 0 * (Laylength[i] - 1)], &temp[2 * 0 * (Laylength[i])], Laylength[i] - 1);//(kx[i+1]-kx[i])^2
		ippsMul_64f_I(&temp[(2 * 0 + 1) * (Laylength[i] - 1)], &temp[(2 * 0 + 1) * (Laylength[i] - 1)], Laylength[i] - 1);//(ky[i+1]-ky[i])^2
		ippsAdd_64f_I(&temp[2 * 0 * (Laylength[i] - 1)], &temp[(2 * 0 + 1) * (Laylength[i] - 1)], Laylength[i] - 1);//(kx[i+1]-kx[i])^2+(ky[i+1]-ky[i])^2
		ippsSqrt_64f_I(&temp[(2 * 0 + 1) * (Laylength[i] - 1)], Laylength[i] - 1);//Sqrt[(kx[i+1]-kx[i])^2+(ky[i+1]-ky[i])^2]
		for (int j = 0; j < Laylength[i]; ++j) {
			temp[j] = 0;

			for (int k = 0; k < j; ++k) {
				temp[j] = temp[j] + temp[(2 * 0 + 1) * (Laylength[i] - 1) + k];//cumulative sum
			}
		}
		for (int j = 0; j < (temp[Laylength[i]-1]/distance); ++j) {
			temp2[j] = distance * j;
		}
		count = 1;
		outx[finalPoint] = kx[cumLaylength];
		outy[finalPoint] = ky[cumLaylength];
		outz[finalPoint] = kz[cumLaylength];
		for (int k = 0; k < Laylength[i]; ++k) {
			if (count > floor(temp[Laylength[i] - 1] / distance)) break;
			if (temp[k] >= temp2[count]) {
				outx[finalPoint + count] = kx[cumLaylength+k - 1] + (kx[cumLaylength+k] - kx[cumLaylength+k - 1])*(temp2[count] - temp[k - 1]) / (temp[k] - temp[k - 1]);
				outz[finalPoint + count] = kz[cumLaylength];
				outy[finalPoint + count] = ky[cumLaylength+k - 1] + (ky[cumLaylength+k] - ky[cumLaylength+k - 1])*(temp2[count] - temp[k - 1]) / (temp[k] - temp[k - 1]);
				++count;

			}
		}
		finalPoint = finalPoint + count;
		cumLaylength = cumLaylength + Laylength[i];
		//cout << finalPoint << endl;
		

	}
	return finalPoint - 1;
}
FindFermi::FindFermi( Ipp64f * param)
{
	int current=0;
	//DataExtractor extractor(name);
	//Ipp64f params[9] = { 0.074, 475, 525, -60, 16, 1000, 0.5, 17, 8 };
	UpdatePar(param);
	fineN = 400;//innitial grid inplane
	NumLeng = 0;//number of contour 

	cdevs =4;//Kz grid
	//Ipp64f * starts = extractor.getDataArray();
	//int nPoints = floor((extractor.getNumberOfLines()) / 2);
	nPoints = 0;
	nfinepoint = (fineN*fineN)*cdevs;
	tol = 1E-10;
	subMaxR = new Ipp64f[1];
	temp1 = new Ipp64f[20 * nfinepoint];
	argCos = new Ipp64f[nfinepoint];
	argSin = new Ipp64f[nfinepoint];
	argkz = new Ipp64f[nfinepoint];
	zeros = new Ipp64f[nfinepoint];
	ippsSet_64f(0,zeros,nfinepoint);
	kz = new Ipp64f[nfinepoint];
	funcval = new Ipp64f[4*nfinepoint];
	tfunc3 = new Ipp64f[nfinepoint];
	tfunc1 = new Ipp64f[nfinepoint];
	tfunc2 = new Ipp64f[nfinepoint];
	tfunc = new Ipp64f[nfinepoint];
	
	lengArr = new int[20*cdevs];
	kx = new Ipp64f[nfinepoint];
	ky = new Ipp64f[nfinepoint];
	finDis = 3.1415926 / 3.747665940 / 15;
	tempx1 = new Ipp64f[nfinepoint];
	tempy1 = new Ipp64f[nfinepoint];
	tempx2 = new Ipp64f[nfinepoint];
	tempy2 = new Ipp64f[nfinepoint];
	tempx3 = new Ipp64f[nfinepoint];
	tempy3 = new Ipp64f[nfinepoint];
	tempx4 = new Ipp64f[nfinepoint];
	tempy4 = new Ipp64f[nfinepoint];
	tempz = new Ipp64f[nfinepoint];
	//Ipp64f *startpoint = new Ipp64f[3 * nPoints];
	for (int i = 0; i < cdevs; i++)
	{
		for (int j = 0; j < fineN; j++)
		{
			for(int k = 0;k<fineN;k++)
			{
				tempx1[i*fineN*fineN + k * fineN + j] = -3.1415926 / 3.747665940 + 2 * 3.1415926 / 3.747665940 / fineN  * j;
				tempy1[i*fineN*fineN + k * fineN + j] = -3.1415926 / 3.747665940 + 2 * 3.1415926 / 3.747665940 / fineN * k;
				tempx2[i*fineN*fineN + k * fineN + j] = -3.1415926 / 3.747665940 + 2 * 3.1415926 / 3.747665940 / fineN * (j+1);
				tempy2[i*fineN*fineN + k * fineN + j] = -3.1415926 / 3.747665940 + 2 * 3.1415926 / 3.747665940 / fineN * (k);
				tempx3[i*fineN*fineN + k * fineN + j] = -3.1415926 / 3.747665940 + 2 * 3.1415926 / 3.747665940 / fineN * j;
				tempy3[i*fineN*fineN + k * fineN + j] = -3.1415926 / 3.747665940 + 2 * 3.1415926 / 3.747665940 / fineN * (k+1);
				tempx4[i*fineN*fineN + k * fineN + j] = -3.1415926 / 3.747665940 + 2 * 3.1415926 / 3.747665940 / fineN * (j+1);
				tempy4[i*fineN*fineN + k * fineN + j] = -3.1415926 / 3.747665940 + 2 * 3.1415926 / 3.747665940 / fineN * (k+1);
				tempz[i*fineN*fineN+k*fineN + j] = -2 * 3.1415926 / 13.2 + 4 * 3.1415926 / 13.2 / cdevs / 2 + 4 * 3.1415926 / 13.2 / cdevs * i;
			}
			/*cout << starts[2 * i] << " " << starts[2 * i + 1] << endl;*/
			
			//kz[i*fineN + j] = -2 * 3.1415926 / 13.2  + 4 * 3.1415926 / 13.2 / cdevs * i;
			

		}
	}

	ippsMulC_64f(tempz, 13.2, argkz, nfinepoint);
	cout << "parameter: ";
	for (int i = 0; i < 10; ++i) {
		cout << params[i] << ", ";
	}
	cout << endl;
	func(params, argkz, tempx1, tempy1, nfinepoint, temp1, tfunc);
	func(params, argkz, tempx2, tempy2, nfinepoint, temp1, tfunc1);
	func(params, argkz, tempx3, tempy3, nfinepoint, temp1, tfunc2);
	func(params, argkz, tempx4, tempy4, nfinepoint, temp1, tfunc3);

	ippsMaxEvery_64f(zeros, tfunc, funcval, nfinepoint);
	ippsMaxEvery_64f(zeros, tfunc1, &funcval[nfinepoint], nfinepoint);
	ippsMaxEvery_64f(zeros, tfunc2, &funcval[2 * nfinepoint], nfinepoint);
	ippsMaxEvery_64f(zeros, tfunc3, &funcval[3 * nfinepoint], nfinepoint);

	ippsDiv_64f_I(tfunc, funcval, nfinepoint);
	ippsDiv_64f_I(tfunc1, &funcval[nfinepoint], nfinepoint);
	ippsDiv_64f_I(tfunc2, &funcval[2 * nfinepoint], nfinepoint);
	ippsDiv_64f_I(tfunc3, &funcval[3 * nfinepoint], nfinepoint);

	ippsAdd_64f_I(&funcval[nfinepoint], funcval, nfinepoint);
	ippsAdd_64f_I(&funcval[2 * nfinepoint], funcval, nfinepoint);
	ippsAdd_64f_I(&funcval[3 * nfinepoint], funcval, nfinepoint);
	nPoints = 0;
	NumLeng = 0;
	for (int k = 0; k < cdevs; ++k) {
		for (int j = 0; j < fineN / 2; ++j) {
			for (int i = 0 ; i < fineN*fineN; ++i) {

				if (funcval[i+ k * (fineN*fineN)] < 3.5 && funcval[i+k * (fineN*fineN)] > 0.5 && ( (i % fineN == j) || ((i ) % fineN == (fineN - j -1))|| ((i) / fineN == j)||((i ) / fineN == (fineN - j - 1)))){
					kx[nPoints] = 0.25*(tempx1[i + k * (fineN*fineN)] + tempx2[i + k * (fineN*fineN)] + tempx3[i + k * (fineN*fineN)] + tempx4[i + k * (fineN*fineN)]);
					ky[nPoints] = 0.25*(tempy1[i + k * (fineN*fineN)] + tempy2[i + k * (fineN*fineN)] + tempy3[i + k * (fineN*fineN)] + tempy4[i + k * (fineN*fineN)]);
					kz[nPoints] = tempz[i + k * (fineN*fineN)];
					funcval[i + k * (fineN*fineN)] = 0;
					nPoints = nPoints + 1;
					current = i;
					NumLeng = NumLeng + 1;
					lengArr[NumLeng - 1] = 1;
					cout << lengArr[NumLeng - 1] << endl;
					while (true) {
						if ((current ) % fineN != 0 && funcval[current-1 + k * (fineN*fineN)] < 3.5 && funcval[current-1 + k * (fineN*fineN)] > 0.5) {
							kx[nPoints] = 0.25*(tempx1[current - 1 + k * (fineN*fineN)] + tempx2[current - 1 + k * (fineN*fineN)] + tempx3[current - 1 + k * (fineN*fineN)] + tempx4[current - 1 + k * (fineN*fineN)]);
							ky[nPoints] = 0.25*(tempy1[current - 1 + k * (fineN*fineN)] + tempy2[current - 1 + k * (fineN*fineN)] + tempy3[current - 1 + k * (fineN*fineN)] + tempy4[current - 1 + k * (fineN*fineN)]);
							kz[nPoints] = tempz[current - 1 + k * (fineN*fineN)];
							funcval[current - 1 + k * (fineN*fineN)] = 0;
							nPoints = nPoints + 1;
							current = current-1;
							lengArr[NumLeng - 1]= lengArr[NumLeng - 1]+1;
							continue;
						}
					    if ((current ) % fineN != (fineN - 1) && funcval[current + 1 + k * (fineN*fineN)] < 3.5 && funcval[current + 1 + k * (fineN*fineN)] > 0.5) {
							kx[nPoints] = 0.25*(tempx1[current + 1 + k * (fineN*fineN)] + tempx2[current + 1 + k * (fineN*fineN)] + tempx3[current + 1 + k * (fineN*fineN)] + tempx4[current + 1 + k * (fineN*fineN)]);
							ky[nPoints] = 0.25*(tempy1[current + 1 + k * (fineN*fineN)] + tempy2[current + 1 + k * (fineN*fineN)] + tempy3[current + 1 + k * (fineN*fineN)] + tempy4[current + 1 + k * (fineN*fineN)]);
							kz[nPoints] = tempz[current + 1 + k * (fineN*fineN)];
							funcval[current + 1 + k * (fineN*fineN)] = 0;
							nPoints = nPoints + 1;
							current = current + 1;
							lengArr[NumLeng - 1] = lengArr[NumLeng - 1] + 1;
							continue;
						}
						 if ((current) / fineN != 0 && funcval[current - fineN + k * (fineN*fineN)] < 3.5 && funcval[current - fineN + k * (fineN*fineN)] > 0.5) {
							kx[nPoints] = 0.25*(tempx1[current - fineN + k * (fineN*fineN)] + tempx2[current - fineN + k * (fineN*fineN)] + tempx3[current - fineN + k * (fineN*fineN)] + tempx4[current - fineN + k * (fineN*fineN)]);
							ky[nPoints] = 0.25*(tempy1[current - fineN + k * (fineN*fineN)] + tempy2[current - fineN + k * (fineN*fineN)] + tempy3[current - fineN + k * (fineN*fineN)] + tempy4[current - fineN + k * (fineN*fineN)]);
							kz[nPoints] = tempz[current - fineN + k * (fineN*fineN)];
							funcval[current - fineN + k * (fineN*fineN)] = 0;
							nPoints = nPoints + 1;
							current = current - fineN;
							lengArr[NumLeng - 1] = lengArr[NumLeng - 1] + 1;
							continue;
						}
						if ((current) / fineN != (fineN - 1) && funcval[current + fineN + k * (fineN*fineN)] < 3.5 && funcval[current + fineN + k * (fineN*fineN)] > 0.5) {
							kx[nPoints] = 0.25*(tempx1[current + fineN + k * (fineN*fineN)] + tempx2[current + fineN + k * (fineN*fineN)] + tempx3[current + fineN + k * (fineN*fineN)] + tempx4[current + fineN + k * (fineN*fineN)]);
							ky[nPoints] = 0.25*(tempy1[current + fineN + k * (fineN*fineN)] + tempy2[current + fineN + k * (fineN*fineN)] + tempy3[current + fineN + k * (fineN*fineN)] + tempy4[current + fineN + k * (fineN*fineN)]);
							kz[nPoints] = tempz[current + fineN + k * (fineN*fineN)];
							funcval[current + fineN + k * (fineN*fineN)] = 0;
							nPoints = nPoints + 1;
							current = current + fineN;
							lengArr[NumLeng - 1] = lengArr[NumLeng - 1] + 1;
							continue;
						}
						//four corners
						if ((current) % fineN != 0 && (current) / fineN != 0 && funcval[current - 1 - fineN + k * (fineN*fineN)] < 3.5 && funcval[current - 1 - fineN + k * (fineN*fineN)] > 0.5) {
							kx[nPoints] = 0.25*(tempx1[current - 1 - fineN + k * (fineN*fineN)] + tempx2[current - 1 - fineN + k * (fineN*fineN)] + tempx3[current - 1 - fineN + k * (fineN*fineN)] + tempx4[current - 1 - fineN + k * (fineN*fineN)]);
							ky[nPoints] = 0.25*(tempy1[current - 1 - fineN + k * (fineN*fineN)] + tempy2[current - 1 - fineN + k * (fineN*fineN)] + tempy3[current - 1 - fineN + k * (fineN*fineN)] + tempy4[current - 1 - fineN + k * (fineN*fineN)]);
							kz[nPoints] = tempz[current - 1 - fineN + k * (fineN*fineN)];
							funcval[current - 1 - fineN + k * (fineN*fineN)] = 0;
							nPoints = nPoints + 1;
							current = current - 1 - fineN;
							lengArr[NumLeng - 1] = lengArr[NumLeng - 1] + 1;
							continue;
						}
						if ((current) % fineN != (fineN - 1) && (current) / fineN != 0 &&funcval[current + 1 - fineN + k * (fineN*fineN)] < 3.5 && funcval[current + 1 - fineN + k * (fineN*fineN)] > 0.5) {
							kx[nPoints] = 0.25*(tempx1[current + 1 - fineN + k * (fineN*fineN)] + tempx2[current + 1 - fineN + k * (fineN*fineN)] + tempx3[current + 1 - fineN + k * (fineN*fineN)] + tempx4[current + 1 - fineN + k * (fineN*fineN)]);
							ky[nPoints] = 0.25*(tempy1[current + 1 - fineN + k * (fineN*fineN)] + tempy2[current + 1 - fineN + k * (fineN*fineN)] + tempy3[current + 1 - fineN + k * (fineN*fineN)] + tempy4[current + 1 - fineN + k * (fineN*fineN)]);
							kz[nPoints] = tempz[current + 1 - fineN + k * (fineN*fineN)];
							funcval[current + 1 - fineN + k * (fineN*fineN)] = 0;
							nPoints = nPoints + 1;
							current = current + 1 - fineN;
							lengArr[NumLeng - 1] = lengArr[NumLeng - 1] + 1;
							continue;
						}
						if ((current) / fineN != (fineN - 1)&& (current) % fineN != 0 && funcval[current - 1 + fineN + k * (fineN*fineN)] < 3.5 && funcval[current + fineN - 1 + k * (fineN*fineN)] > 0.5) {
							kx[nPoints] = 0.25*(tempx1[current + fineN - 1 + k * (fineN*fineN)] + tempx2[current + fineN - 1 + k * (fineN*fineN)] + tempx3[current + fineN - 1 + k * (fineN*fineN)] + tempx4[current + fineN - 1 + k * (fineN*fineN)]);
							ky[nPoints] = 0.25*(tempy1[current + fineN - 1 + k * (fineN*fineN)] + tempy2[current + fineN - 1 + k * (fineN*fineN)] + tempy3[current + fineN - 1 + k * (fineN*fineN)] + tempy4[current + fineN - 1 + k * (fineN*fineN)]);
							kz[nPoints] = tempz[current + fineN - 1 + k * (fineN*fineN)];
							funcval[current + fineN - 1 + k * (fineN*fineN)] = 0;
							nPoints = nPoints + 1;
							current = current + fineN - 1;
							lengArr[NumLeng - 1] = lengArr[NumLeng - 1] + 1;
							continue;
						}
						if ((current) / fineN != (fineN - 1) && (current) % fineN != (fineN - 1) && funcval[current + fineN + 1 + k * (fineN*fineN)] < 3.5 && funcval[current + fineN + 1 + k * (fineN*fineN)] > 0.5) {
							kx[nPoints] = 0.25*(tempx1[current + fineN + 1 + k * (fineN*fineN)] + tempx2[current + fineN + 1 + k * (fineN*fineN)] + tempx3[current + fineN + 1 + k * (fineN*fineN)] + tempx4[current + fineN + 1 + k * (fineN*fineN)]);
							ky[nPoints] = 0.25*(tempy1[current + fineN + 1 + k * (fineN*fineN)] + tempy2[current + fineN + 1 + k * (fineN*fineN)] + tempy3[current + fineN + 1 + k * (fineN*fineN)] + tempy4[current + fineN + 1 + k * (fineN*fineN)]);
							kz[nPoints] = tempz[current + fineN + 1 + k * (fineN*fineN)];
							funcval[current + fineN + 1 + k * (fineN*fineN)] = 0;
							nPoints = nPoints + 1;
							current = current + fineN + 1;
							lengArr[NumLeng - 1] = lengArr[NumLeng - 1] + 1;
							continue;
						}
						cout << lengArr[NumLeng - 1] << endl;
						break;

					}
					//cout << ((i - k * (fineN*fineN)) / fineN == (fineN - j - 1) )<< endl;

				}
			}
		}
	}
	nPoints = interfunc(kx, ky, kz, lengArr, NumLeng, finDis, temp1, tempx1, tempx2, tempy2, tempz);
	//cout << nPoints << "Second times	" << NumLeng <<"	"<< endl;
}

int FindFermi::UpdatePar(Ipp64f * param)
{
	for (int i = 0; i < 10; ++i) {
		params[i] = param[i];
	}

	return 0;
}

int FindFermi::PrintPar()
{
	cout << "parameter: ";
	for (int i = 0; i < 10; ++i) {
		cout << params[i] << ", ";
	}
	cout << endl;
	return 0;
}



int FindFermi::ReturnStart(Ipp64f * startpoint)
{
	
		




		
		


	
	//cout << nPoints << "Second times	" << NumLeng << "	" << endl;





	for (int i = 0; i < nPoints; i++) {
		startpoint[i * 3] = tempx2[i];
		startpoint[i * 3 + 1] = tempy2[i];
		startpoint[i * 3 + 2] = tempz[i];
	}
	ofstream fout;
	fout.open("FindFermi.dat");
	fout.precision(15);

	for (int i = 0; i < nPoints; i++) {

		fout << startpoint[i * 3] << "\t" << startpoint[i * 3 + 1] << "\t" << startpoint[i * 3 + 2] << endl;
	//	fout << theta[i] << "\t" << r[i] << "\t" << kz[i] << endl;
		//cout << kx[i] << "\t" << ky[i] << "\t" << kz[i] << endl;
	}

	fout.close();
	//while (true);
	return 0;
}

int FindFermi::ReturnNumPoint()
{
	cout << "nPoints = " << &nPoints << "  " << nPoints << endl;
	return nPoints;
}


FindFermi::~FindFermi()
{
	delete subMaxR;
	delete temp1;
	delete argCos;
	delete argSin;
	delete argkz;
	delete kz;
	delete funcval;
	delete tfunc1;
	delete tfunc2;
	delete tfunc;
	delete tfunc3;
	delete tempx1;
	delete tempy1;
	delete tempx2;
	delete tempy2;
	delete tempx3;
	delete tempy3;
	delete tempx4;
	delete tempy4;
	delete tempz;
	delete lengArr;
	delete kx;
	delete ky;

}
