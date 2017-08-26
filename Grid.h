#ifndef GRID_H
#define GRID_H

#include <iostream> // and I put this here because.......
#include <random>

#include "CTextures.h"
#include "SpaceCamera.h"
#include "hgtreader.h"

//#include "Vec3.h" // 3D vector class

// Currently replacing the old Normal3d and Point3d classes with a combined Vec3 class.
// This should hopefully cut down on code duplication.
// Undecided on whether to have another class which then combines these two classes 
// or to simple create two instances in the Grid.h class.

/* 
class VertexData() {
    Vec3 point3d;
    Vec3 normal3d;
    
    or 
    
    std::vector<Vec3> point3d;
    std::vector<Vec3> normal3d;
    
    currently favouring the first option, as it seams to give more flexability.
    addtovector(Vec3& v) {
        v.push_back(Vec3());
    }
}
*/

// TO DO: 22/08/14
/*
   As the old point3d and normal3d classes are being superseded by the new Vec3 class 
   a lot of the normal calculation code will be made superflicious. 
   The intention is to replace this with a new implemntation which uses vector math
   functions in the Vec3 class. This should make things a little less haphazed.
*/


class Normal3d {
private:
	double x;
	double y;
	double z;
public:
	Normal3d(double nx, double ny, double nz) : x(nx), y(ny), z(nz) {
	}
	Normal3d() {};
	~Normal3d() {}
	double getx()  { return x; }
    double gety()  { return y; }
    double getz()  { return z; }
};

class Point3d {
private:
    double x;
    double y;
    double z;

public:
	Point3d(double nx, double ny, double nz) : x(nx), y(ny), z(nz) {		
	}
	Point3d() {};
	~Point3d() {}

	double getx()  { return x; }
    double gety()  { return y; }
    double getz()  { return z; }

	Point3d operator+(Point3d a) {
		Point3d b;

		b.x = this->x + a.x;
		b.y = this->y + a.y;
		b.z = this->z + a.z;
		return b;
	}

	Point3d operator-(Point3d a) {
		Point3d b;

		b.x = this->x - a.x;
		b.y = this->y - a.y;
		b.z = this->z - a.z;
		return b;
	}

	Point3d operator/(double a) {
		Point3d b;

		b.x = this->x / a;
		b.y = this->y / a;
		b.z = this->z / a;
		return b;
	}	

	Normal3d normal[1];
	
};

class RandomNumbers {
private:
	std::random_device generator;
public:
	RandomNumbers() {
		
	}

	double GenerateUniformRealDistribution_Range(double displacement_factor, double input) {

		double range = input * displacement_factor;

		double min = 0 - range;
		double max = 0 + range;

		std::uniform_real_distribution<> distribution{ min, max };

		return distribution(generator);
	}

	double GenerateUniformRealDistribution(double min, double max) {

		std::uniform_real_distribution<> distribution{ min, max };

		return distribution(generator);
	}

	double GenerateNormalDistribution() {

		std::normal_distribution<double> distribution(0.0, 0.5);
		return distribution(generator);

	}

	double GenerateChiSquaredDistribution(double num) {

		std::chi_squared_distribution<double> distribution(num);
		return distribution(generator);

	}

	int GenerateUniformIntDistribution(int min, int max) {

		std::uniform_int_distribution<> distribution{ min, max };

		return distribution(generator);
	}

	Point3d GenerateForVector() {
		double gen_1 = GenerateUniformRealDistribution(-0.025, 0.025);
		Point3d vec_1(0, 0, 0);

		double gen_2 = GenerateNormalDistribution();
		Point3d vec_2(0, 0, 0);

		double gen_3 = GenerateChiSquaredDistribution(0.025);
		Point3d vec_3(0, 0, 0);

		Point3d vec_4 = (vec_1 + vec_3) / 2.0;

		return vec_4;
	}
};

class NoiseContainer {
private:
	double frequancy;
	double amplitude;

	int       length;
	int       rangevalues;

	std::vector<double> points;

	std::vector<double> xvalues;
	std::vector<double> yvalues;
	std::vector<double> zvalues;

public:
	NoiseContainer(double nfrequancy, double namplitude, int nlength, int nrangevalues) :

		frequancy(nfrequancy),
		amplitude(namplitude),
		length(nlength),
		rangevalues(nrangevalues)
	{

		NoiseFunc();

	}

	~NoiseContainer() {}

	void NoiseFunc() {

		RandomNumbers rn;

		double initialvalue_x = rn.GenerateUniformRealDistribution(0.0, 2.0 * PI);
		double initialvalue_y = rn.GenerateUniformRealDistribution(0.0, 2.0 * PI);

		xvalues.push_back(initialvalue_x);
		yvalues.push_back(initialvalue_y);

		double variance = 0.05; // 0.025
		
		switch (rangevalues) {

			case 0:
				for (int i = 1; i < length; i++) {
					double newvalue_x = xvalues[i - 1] + rn.GenerateUniformRealDistribution(0.0, variance);
					xvalues.push_back(newvalue_x);

					double newvalue_y = yvalues[i - 1] + rn.GenerateUniformRealDistribution(0.0, variance);
					yvalues.push_back(newvalue_y);
				}
				break;

			case 1:
				for (int i = 1; i < length; i++) {
					
					double newvalue_x = xvalues[i - 1] + rn.GenerateUniformRealDistribution(-variance, 0.0);
					xvalues.push_back(newvalue_x);
					
					double newvalue_y = yvalues[i - 1] + rn.GenerateUniformRealDistribution(0.0, variance);
					yvalues.push_back(newvalue_y);
				}
				break;

			case 2:
				for (int i = 1; i < length; i++) {
					
					double newvalue_x = xvalues[i - 1] + rn.GenerateUniformRealDistribution(-variance, 0.0);
					xvalues.push_back(newvalue_x);
										
					double newvalue_y = yvalues[i - 1] + rn.GenerateUniformRealDistribution(-variance, 0.0);
					yvalues.push_back(newvalue_y);
				}
				break;

			case 3:
				for (int i = 1; i < length; i++) {
					
					double newvalue_x = xvalues[i - 1] + rn.GenerateUniformRealDistribution(0.0, variance);
					xvalues.push_back(newvalue_x);
										
					double newvalue_y = yvalues[i - 1] + rn.GenerateUniformRealDistribution(-variance, 0.0);
					yvalues.push_back(newvalue_y);
				}
				break;

			default:
				break;

		}


		for (int h = 0; h < yvalues.size(); h++) {

			for (int w = 0; w < xvalues.size(); w++) {

				double newvalue_z = (getYValue(h) * getXValue(w)) / (2.0 * PI);				

				zvalues.push_back(newvalue_z);

			}

		}

		double offset = rn.GenerateUniformIntDistribution(0, 6);

		for (int i = 0; i < zvalues.size(); i++) {
			
			double  noise = amplitude *  sin((getZValue(i) * frequancy) + offset);
			
			points.push_back(noise);

		}

	}

	double getNoiseValue(int index)    { return points[index]; }

	double getXValue(int index)   { return xvalues[index]; }
	double getYValue(int index)   { return yvalues[index]; }
	double getZValue(int index)   { return zvalues[index]; }

	double getLength()            { return length; }

};

class Triangle3d {
private:
	int v1;
    int v2;
    int v3;
	
public:
	Triangle3d(int nv1, int nv2, int nv3) : v1(nv1), v2(nv2), v3(nv3) {
		
	}
	~Triangle3d() {}        
	
	int getv1()  { return v1; }
    int getv2()  { return v2; }
    int getv3()  { return v3; }	
	
};

class CGrid {
private:
	int length;

	int wc;
	int heightcounter;

	

	RandomNumbers rn;
	double frequancy;
	double amplitude;
		

public:
	CGrid(int nlength) : length(nlength) {

		wc = 0;
		heightcounter = 0;
		frequancy = 1.0;
		amplitude = 1.0;		

		initializegrid();
		
	}

	~CGrid() {}

	std::vector<Point3d> point3d;
	std::vector<Triangle3d> triangle3d;	

	std::vector<double> waterheight; // should probably be part of a class that combines all data about a spot on the map	
	std::vector<double> sedimentheight; 

	std::vector<NoiseContainer> noisecontainer;
	std::vector<double> sumofnoise;

	std::vector<int> tempids;
	std::vector<Point3d> tempnormals;		

	std::vector<double> noise;

	

	void addPoint3d(double nx, double ny, double nz) { point3d.push_back(Point3d(nx, ny, nz)); }
	void addTriangle3d(int nv1, int nv2, int nv3)    { triangle3d.push_back(Triangle3d(nv1, nv2, nv3)); }

	void initializegrid() {
		
		for (int i = 0; i < length * length; i++) {
			sumofnoise.push_back(0.0);
			//waterheight.push_back(9.600); // 0.004
			//sedimentheight.push_back(0.0);
		}

		double max = 0;
		double min = 0;

		for (int i = 0; i < 4; i++) {			
			noisecontainer.push_back(NoiseContainer(frequancy, amplitude, length, 0));
			noisecontainer.push_back(NoiseContainer(frequancy, amplitude, length, 1));
			noisecontainer.push_back(NoiseContainer(frequancy, amplitude, length, 2));
			noisecontainer.push_back(NoiseContainer(frequancy, amplitude, length, 3));					

			for (int i_point = 0; i_point < length * length; i_point++) {

				double tempsum = 0;

				for (int i_container = noisecontainer.size() - 4; i_container < noisecontainer.size(); i_container++) {
					tempsum = tempsum + noisecontainer[i_container].getNoiseValue(i_point);	
				}

				tempsum = tempsum / 4.0;					
				sumofnoise[i_point] = sumofnoise[i_point] + tempsum;
			}

			frequancy = frequancy * 2.0;
			amplitude = amplitude / 2.0;			
		}			

		for (int h = 0; h < length; h++) {
			for (int w = 0; w < length; w++) {
				addPoint3d(w, 0.0, h);			
			}
		}				
		
		

		// cosine interpolation
		for (int p = 0; p < 6; p++) {
			int hc = 1;
			for (int i = 1; i < sumofnoise.size() - 1; i++) {
				double mu = (1 - cos(0.5 * PI)) / 2.0;

				double v1 = sumofnoise[i] * (1 - mu) + sumofnoise[i + 1] * mu;
				double v2 = sumofnoise[i - 1] * (1 - mu) + sumofnoise[i] * mu;

				double v3 = v2 * (1 - mu) + v1 * mu;
				sumofnoise[i] = v3;

				if (i == (length * hc) - 2) {
					i = i + 2;
					hc++;
				}
			}
			for (int i = length; i < sumofnoise.size() - length; i++) {
				double mu = (1 - cos(0.5 * PI)) / 2.0;

				double v1 = (sumofnoise[i] * (1 - mu) + sumofnoise[i + length] * mu);
				double v2 = (sumofnoise[i - length] * (1 - mu) + sumofnoise[i] * mu);

				double v3 = v2 * (1 - mu) + v1 * mu;
				sumofnoise[i] = v3;
			}
		}

		// normalise the noise to range 0 - 1
		for (int i = 0; i < sumofnoise[i]; i++) {
			if (sumofnoise[i] > max) { max = sumofnoise[i]; }
			if (sumofnoise[i] < min) { min = sumofnoise[i]; }
		}
		for (int i = 0; i < sumofnoise[i]; i++) {
			sumofnoise[i] = (sumofnoise[i] - min) / (max - min);
		}
		// assign points and scale the mesh
		for (int i = 0; i < length * length; i++) {			
			point3d[i] = Point3d(point3d[i].getx(), sumofnoise[i], point3d[i].getz());
			point3d[i] = Point3d(point3d[i].getx(), 50 * pow(point3d[i].gety(), 2.0), point3d[i].getz());
			//waterheight[i] = waterheight[i]; //set to the same scale as the terrain
		}

		/*for (int i = 0; i < 4000; i++) {
			calculatewatermovement();			
		}*/

		/*for (int i = 0; i < waterheight.size(); i++) {
			waterheight[i] = 0.00125;			
		}*/
		//calculatewatermovement();
		waterheight.clear();
		sedimentheight.clear();

		refinetriangles(); 		
		vertexnormals();
	}

	void refinetriangles() {

		int trinumber = ((length - 1) * (length - 1)); // number of triangles
		int index = 0;
		int h = 1;

		for (int a = 0; a < trinumber; a++) {

			//addTriangle3d(index, index + 1, index + gridwidth+1);
			//addTriangle3d(index + gridwidth + 1, index + gridwidth, index);

			addTriangle3d(index + length + 1, index + 1, index);
			addTriangle3d(index, index + length, index + length + 1);

			if (index == (length * h) - 2) {
				index = index + 2;
				h++;
			}
			else { index = index + 1; }
		}
	}

	void findlowestneighbour(std::vector<double> v, double &nmax, int &npos) { 

		// find the greatest differance
		for (int i = 0; i < v.size(); i++) {
			if (v[i] > 0.0) {
				if (v[i] > nmax) {
					nmax = v[i];
					// record its position
					npos = i;
				}
			}
			// else the neighbour is higher than starting cell.
			else { }
		}
	}

	void transportwater(double &na, double &nb, double transferrate, double erosionrate, int pos, int wc, int wc2) {

		if (na >= transferrate) {
			if (point3d[wc].gety() - erosionrate >= point3d[wc2].gety()) {
				point3d[wc] = Point3d(point3d[wc].getx(), point3d[wc].gety() - erosionrate, point3d[wc].getz());
				sedimentheight[wc2] = sedimentheight[wc2] + erosionrate;
			}
			else {
				// transfer half the differance
				point3d[wc] = Point3d(point3d[wc].getx(), point3d[wc].gety() - (point3d[wc].gety() - point3d[wc2].gety()), point3d[wc].getz());
				sedimentheight[wc2] = sedimentheight[wc2] + (point3d[wc].gety() - point3d[wc2].gety());
			}
			nb = nb + transferrate;
			na = na - transferrate;			
		}
		else if (na < transferrate) {

			if (na <= 0) { na = 0; }
			if (na > 0.0) {
				// find percentage of errosian rate
				double adjustederrosionrate = (na / transferrate) * erosionrate;

				if (point3d[wc].gety() - adjustederrosionrate >= point3d[wc2].gety()) {
					point3d[wc] = Point3d(point3d[wc].getx(), point3d[wc].gety() - adjustederrosionrate, point3d[wc].getz());
					sedimentheight[wc2] = sedimentheight[wc2] + adjustederrosionrate;
				}
				else {
					point3d[wc] = Point3d(point3d[wc].getx(), point3d[wc].gety() - (point3d[wc].gety() - point3d[wc2].gety()), point3d[wc].getz());
					sedimentheight[wc2] = sedimentheight[wc2] + (point3d[wc].gety() - point3d[wc2].gety());
				}
			}
			nb = nb + na;
			na = 0.0;			
		}			

	}
	void calculatewatermovement() {
		// wc relpaced a, which was a local in var
		// h used to be another loop like w

		int gw = length - 1;		
		int h = heightcounter;

		std::vector<double> combinedheight;

		double max = 0;
		int pos = -1;

		double erosionrate       = 0.1250;    // 0.00125;
		double sedimentratio     = 0.0001;       // 0.05
		double watertransferrate = 0.001250;      // 0.25
		double evapourationrate  = 0.000312500;  // 0.0031250;

		for (int w = 0; w < length; w++) {

			if (waterheight[wc] == 0) {
				// do nothing 						
			}
			else {
				if ((h == 0) && (w == 0)) {

					combinedheight.push_back((point3d[wc].gety() + waterheight[wc]) -
						(point3d[wc + length].gety() + waterheight[wc + length]));

					combinedheight.push_back((point3d[wc].gety() + waterheight[wc]) -
						(point3d[wc + 1].gety() + waterheight[wc + 1]));

					findlowestneighbour(combinedheight, max, pos);
					// if pos == - 1 then no erosion takes place
					// else home tile gets eroded

					// errode a as sediment and transfer it in water to lowestneighbour
					if (pos == 0) {
						transportwater(waterheight[wc], waterheight[wc + length], watertransferrate, erosionrate, pos, wc, wc + length);
					}
					else if (pos == 1) {
						transportwater(waterheight[wc], waterheight[wc + 1], watertransferrate, erosionrate, pos, wc, wc + 1);
					}

					// Deposit sediment
					if (sedimentheight[wc] > (waterheight[wc] * sedimentratio)) {
						double nonsoluablesediment = sedimentheight[wc] - (waterheight[wc] * sedimentratio);
						point3d[wc] = Point3d(point3d[wc].getx(), point3d[wc].gety() + nonsoluablesediment, point3d[wc].getz());
					}

				}

				else if ((h == 0) && (w < gw)) {

					combinedheight.push_back((point3d[wc].gety() + waterheight[wc]) -
						(point3d[wc + length].gety() + waterheight[wc + length]));

					combinedheight.push_back((point3d[wc].gety() + waterheight[wc]) -
						(point3d[wc + 1].gety() + waterheight[wc + 1]));

					combinedheight.push_back((point3d[wc].gety() + waterheight[wc]) -
						(point3d[wc - 1].gety() + waterheight[wc - 1]));

					findlowestneighbour(combinedheight, max, pos);

					if (pos == 0) {
						transportwater(waterheight[wc], waterheight[wc + length], watertransferrate, erosionrate, pos, wc, wc + length);
					}
					else if (pos == 1) {
						transportwater(waterheight[wc], waterheight[wc + 1], watertransferrate, erosionrate, pos, wc, wc + 1);
					}
					else if (pos == 2) {
						transportwater(waterheight[wc], waterheight[wc - 1], watertransferrate, erosionrate, pos, wc, wc - 1);
					}

					// Deposit sediment
					if (sedimentheight[wc] > (waterheight[wc] * sedimentratio)) {
						double nonsoluablesediment = sedimentheight[wc] - (waterheight[wc] * sedimentratio);
						point3d[wc] = Point3d(point3d[wc].getx(), point3d[wc].gety() + nonsoluablesediment, point3d[wc].getz());
					}

				}
				else if ((h == 0) && (w = gw)) {

					combinedheight.push_back((point3d[wc].gety() + waterheight[wc]) -
						(point3d[wc + length].gety() + waterheight[wc + length]));

					combinedheight.push_back((point3d[wc].gety() + waterheight[wc]) -
						(point3d[wc - 1].gety() + waterheight[wc - 1]));

					findlowestneighbour(combinedheight, max, pos);

					if (pos == 0) {
						transportwater(waterheight[wc], waterheight[wc + length], watertransferrate, erosionrate, pos, wc, wc + length);
					}
					else if (pos == 1) {
						transportwater(waterheight[wc], waterheight[wc - 1], watertransferrate, erosionrate, pos, wc, wc - 1);
					}

					// Deposit sediment
					if (sedimentheight[wc] > (waterheight[wc] * sedimentratio)) {
						double nonsoluablesediment = sedimentheight[wc] - (waterheight[wc] * sedimentratio);
						point3d[wc] = Point3d(point3d[wc].getx(), point3d[wc].gety() + nonsoluablesediment, point3d[wc].getz());
					}
				}
				else if ((h > 0) && (h < gw) && (w == 0)) {

					combinedheight.push_back((point3d[wc].gety() + waterheight[wc]) -
						(point3d[wc + length].gety() + waterheight[wc + length]));

					combinedheight.push_back((point3d[wc].gety() + waterheight[wc]) -
						(point3d[wc + 1].gety() + waterheight[wc + 1]));

					combinedheight.push_back((point3d[wc].gety() + waterheight[wc]) -
						(point3d[wc - length].gety() + waterheight[wc - length]));

					findlowestneighbour(combinedheight, max, pos);

					if (pos == 0) {
						transportwater(waterheight[wc], waterheight[wc + length], watertransferrate, erosionrate, pos, wc, wc + length);
					}
					else if (pos == 1) {
						transportwater(waterheight[wc], waterheight[wc + 1], watertransferrate, erosionrate, pos, wc, wc + 1);
					}
					else if (pos == 2) {
						transportwater(waterheight[wc], waterheight[wc - length], watertransferrate, erosionrate, pos, wc, wc - length);
					}

					// Deposit sediment
					if (sedimentheight[wc] > (waterheight[wc] * sedimentratio)) {
						double nonsoluablesediment = sedimentheight[wc] - (waterheight[wc] * sedimentratio);
						point3d[wc] = Point3d(point3d[wc].getx(), point3d[wc].gety() + nonsoluablesediment, point3d[wc].getz());
					}
				}
				else if ((h > 0) && (h < gw) && (w < gw)) {

					combinedheight.push_back((point3d[wc].gety() + waterheight[wc]) -
						(point3d[wc + length].gety() + waterheight[wc + length]));

					combinedheight.push_back((point3d[wc].gety() + waterheight[wc]) -
						(point3d[wc + 1].gety() + waterheight[wc + 1]));

					combinedheight.push_back((point3d[wc].gety() + waterheight[wc]) -
						(point3d[wc - length].gety() + waterheight[wc - length]));

					combinedheight.push_back((point3d[wc].gety() + waterheight[wc]) -
						(point3d[wc - 1].gety() + waterheight[wc - 1]));

					findlowestneighbour(combinedheight, max, pos);

					if (pos == 0) {
						transportwater(waterheight[wc], waterheight[wc + length], watertransferrate, erosionrate, pos, wc, wc + length);
					}
					else if (pos == 1) {
						transportwater(waterheight[wc], waterheight[wc + 1], watertransferrate, erosionrate, pos, wc, wc + 1);
					}
					else if (pos == 2) {
						transportwater(waterheight[wc], waterheight[wc - length], watertransferrate, erosionrate, pos, wc, wc - length);
					}
					else if (pos == 3) {
						transportwater(waterheight[wc], waterheight[wc - 1], watertransferrate, erosionrate, pos, wc, wc - 1);
					}

					// Deposit sediment
					if (sedimentheight[wc] > (waterheight[wc] * sedimentratio)) {
						double nonsoluablesediment = sedimentheight[wc] - (waterheight[wc] * sedimentratio);
						point3d[wc] = Point3d(point3d[wc].getx(), point3d[wc].gety() + nonsoluablesediment, point3d[wc].getz());
					}
				}
				else if ((h > 0) && (h < gw) && (w == gw)) {

					combinedheight.push_back((point3d[wc].gety() + waterheight[wc]) -
						(point3d[wc + length].gety() + waterheight[wc + length]));

					combinedheight.push_back((point3d[wc].gety() + waterheight[wc]) -
						(point3d[wc - length].gety() + waterheight[wc - length]));

					combinedheight.push_back((point3d[wc].gety() + waterheight[wc]) -
						(point3d[wc - 1].gety() + waterheight[wc - 1]));

					findlowestneighbour(combinedheight, max, pos);

					if (pos == 0) {
						transportwater(waterheight[wc], waterheight[wc + length], watertransferrate, erosionrate, pos, wc, wc + length);
					}
					else if (pos == 1) {
						transportwater(waterheight[wc], waterheight[wc - length], watertransferrate, erosionrate, pos, wc, wc - length);
					}
					else if (pos == 2) {
						transportwater(waterheight[wc], waterheight[wc - 1], watertransferrate, erosionrate, pos, wc, wc - 1);
					}

					// Deposit sediment
					if (sedimentheight[wc] > (waterheight[wc] * sedimentratio)) {
						double nonsoluablesediment = sedimentheight[wc] - (waterheight[wc] * sedimentratio);
						point3d[wc] = Point3d(point3d[wc].getx(), point3d[wc].gety() + nonsoluablesediment, point3d[wc].getz());
					}
				}
				else if ((h == gw) && (w == 0)) {

					combinedheight.push_back((point3d[wc].gety() + waterheight[wc]) -
						(point3d[wc + 1].gety() + waterheight[wc + 1]));

					combinedheight.push_back((point3d[wc].gety() + waterheight[wc]) -
						(point3d[wc - length].gety() + waterheight[wc - length]));

					findlowestneighbour(combinedheight, max, pos);

					if (pos == 0) {
						transportwater(waterheight[wc], waterheight[wc + 1], watertransferrate, erosionrate, pos, wc, wc + 1);
					}
					else if (pos == 1) {
						transportwater(waterheight[wc], waterheight[wc - length], watertransferrate, erosionrate, pos, wc, wc - length);
					}

					// Deposit sediment
					if (sedimentheight[wc] > (waterheight[wc] * sedimentratio)) {
						double nonsoluablesediment = sedimentheight[wc] - (waterheight[wc] * sedimentratio);
						point3d[wc] = Point3d(point3d[wc].getx(), point3d[wc].gety() + nonsoluablesediment, point3d[wc].getz());
					}
				}
				else if ((h == gw) && (w > 0) && (w < gw)) {

					combinedheight.push_back((point3d[wc].gety() + waterheight[wc]) -
						(point3d[wc + 1].gety() + waterheight[wc + 1]));

					combinedheight.push_back((point3d[wc].gety() + waterheight[wc]) -
						(point3d[wc - length].gety() + waterheight[wc - length]));

					combinedheight.push_back((point3d[wc].gety() + waterheight[wc]) -
						(point3d[wc - 1].gety() + waterheight[wc - 1]));

					findlowestneighbour(combinedheight, max, pos);

					if (pos == 0) {
						transportwater(waterheight[wc], waterheight[wc + 1], watertransferrate, erosionrate, pos, wc, wc + 1);
					}
					else if (pos == 1) {
						transportwater(waterheight[wc], waterheight[wc - length], watertransferrate, erosionrate, pos, wc, wc - length);
					}
					else if (pos == 2) {
						transportwater(waterheight[wc], waterheight[wc - 1], watertransferrate, erosionrate, pos, wc, wc - 1);
					}

					// Deposit sediment
					if (sedimentheight[wc] > (waterheight[wc] * sedimentratio)) {
						double nonsoluablesediment = sedimentheight[wc] - (waterheight[wc] * sedimentratio);
						point3d[wc] = Point3d(point3d[wc].getx(), point3d[wc].gety() + nonsoluablesediment, point3d[wc].getz());
					}
				}
				else if ((h == gw) && (w == gw)) {

					combinedheight.push_back((point3d[wc].gety() + waterheight[wc]) -
						(point3d[wc - length].gety() + waterheight[wc - length]));

					combinedheight.push_back((point3d[wc].gety() + waterheight[wc]) -
						(point3d[wc - 1].gety() + waterheight[wc - 1]));

					findlowestneighbour(combinedheight, max, pos);

					if (pos == 0) {
						transportwater(waterheight[wc], waterheight[wc - length], watertransferrate, erosionrate, pos, wc, wc - length);
					}
					else if (pos == 1) {
						transportwater(waterheight[wc], waterheight[wc - 1], watertransferrate, erosionrate, pos, wc, wc - 1);
					}

					// Deposit sediment
					if (sedimentheight[wc] > (waterheight[wc] * sedimentratio)) {
						double nonsoluablesediment = sedimentheight[wc] - (waterheight[wc] * sedimentratio);
						point3d[wc] = Point3d(point3d[wc].getx(), point3d[wc].gety() + nonsoluablesediment, point3d[wc].getz());
					}
				}

				//if (waterheight[wc] < 0.0) { waterheight[wc] = 0.0; }
				//if (sedimentheight[wc] < 0.0) { sedimentheight[wc] = 0.0; }

				combinedheight.clear();
				max = 0;
				pos = -1;
			} // else statment

			//waterheight[wc] = waterheight[wc] - evapourationrate;
			//if (waterheight[wc] < 0.0) { waterheight[wc] = 0.0; }

			combinedheight.clear();
			max = 0;
			pos = -1;			
			wc++;
		} // for loop w	
		
		
		heightcounter++;
		if (heightcounter == length) { 
			heightcounter = 0; 
			wc = 0;
		}
		
	}

	

	void vertexnormals() {
		int gw = length - 1;
		int a = 0;
		for (int h = 0; h < length; h++) {
			for (int w = 0; w < length; w++) {				

				double fx = 0; double fy = 0; double fz = 0;

				if ((h == 0) && (w == 0)) {
					// bottom row far left
					tempids.push_back(a);				// 0
					tempids.push_back(a + gw + 1);   // 1
					tempids.push_back(a + gw + 2);   // 2
					tempids.push_back(a + 1);             // 3

					// Calculate normals for each of the tempids.
					// Specific assignment

					assignnormal(tempids[0], tempids[1], tempids[2]);
					assignnormal(tempids[0], tempids[2], tempids[3]);
					//assignnormal(tempids[0],tempids[3],4);
					//assignnormal(0,4,1);

					// Calculate the Sum of the normals 

					for (int i = 0; i < 2; i++) {
						fx = fx + tempnormals[i].getx();
						fy = fy + tempnormals[i].gety();
						fz = fz + tempnormals[i].getz();
					}

					fx = fx / 2;
					fy = fy / 2;
					fz = fz / 2;

					// Assign final normal to the point.

					point3d[a].normal[0] = Normal3d(fx, fy, fz);					

				}
				else if ((h == 0) && (w < gw)) {
					// bottom row
					tempids.push_back(a);				// 0
					tempids.push_back(a + gw + 1);	// 1
					tempids.push_back(a + gw + 2);	// 2
					tempids.push_back(a + 1);				// 3
					tempids.push_back(a - 1);				// 4
					tempids.push_back(a + gw);		// 5

					// Calculate normals for each of the tempids.
					// Specific assignment

					assignnormal(tempids[0], tempids[1], tempids[2]);
					assignnormal(tempids[0], tempids[2], tempids[3]);
					assignnormal(tempids[0], tempids[4], tempids[5]);
					assignnormal(tempids[0], tempids[5], tempids[1]);

					// Calculate the Sum of the normals 

					for (int i = 0; i < 4; i++) {
						fx = fx + tempnormals[i].getx();
						fy = fy + tempnormals[i].gety();
						fz = fz + tempnormals[i].getz();
					}

					fx = fx / 4;
					fy = fy / 4;
					fz = fz / 4;

					// Assign final normal to the point.
					point3d[a].normal[0] = Normal3d(fx, fy, fz);			
					
				}
				else if ((h == 0) && (w = gw)) {
					// bottom row far right
					tempids.push_back(a);				// 0
					tempids.push_back(a + gw + 1);   // 1
					tempids.push_back(a + gw);		// 2
					tempids.push_back(a - 1);				// 3


					// Calculate normals for each of the tempids.
					// Specific assignment

					assignnormal(tempids[0], tempids[1], tempids[2]);
					assignnormal(tempids[0], tempids[2], tempids[3]);

					// Calculate the Sum of the normals 

					for (int i = 0; i < 2; i++) {
						fx = fx + tempnormals[i].getx();
						fy = fy + tempnormals[i].gety();
						fz = fz + tempnormals[i].getz();
					}

					fx = fx / 2;
					fy = fy / 2;
					fz = fz / 2;

					// Assign final normal to the point.
					point3d[a].normal[0] = Normal3d(fx, fy, fz);
					
				}
				else if ((h>0) && (h < gw) && (w == 0)) {
					// left column
					tempids.push_back(a);				// 0
					tempids.push_back(a + gw + 1);	// 1
					tempids.push_back(a + gw + 2);	// 2
					tempids.push_back(a + 1);				// 3
					tempids.push_back(a - gw);		// 4
					tempids.push_back(a - gw - 1);	// 5

					// Calculate normals for each of the tempids.
					// Specific assignment

					assignnormal(tempids[0], tempids[1], tempids[2]);
					assignnormal(tempids[0], tempids[2], tempids[3]);
					assignnormal(tempids[0], tempids[3], tempids[4]);
					assignnormal(tempids[0], tempids[4], tempids[5]);

					// Calculate the Sum of the normals 

					for (int i = 0; i < 4; i++) {
						fx = fx + tempnormals[i].getx();
						fy = fy + tempnormals[i].gety();
						fz = fz + tempnormals[i].getz();
					}

					fx = fx / 4;
					fy = fy / 4;
					fz = fz / 4;

					// Assign final normal to the point.
					point3d[a].normal[0] = Normal3d(fx, fy, fz);

				}
				else if ((h>0) && (h < gw) && (w < gw)) {
					// everything in the middle
					tempids.push_back(a);				// 0
					tempids.push_back(a + gw + 1);	// 1
					tempids.push_back(a + gw + 2);	// 2
					tempids.push_back(a + 1);				// 3
					tempids.push_back(a - gw);		// 4
					tempids.push_back(a - gw - 1);	// 5
					tempids.push_back(a - gw - 2);	// 6
					tempids.push_back(a - 1);				// 7
					tempids.push_back(a + gw);		// 8

					// Calculate normals for each of the tempids.
					// Specific assignment

					assignnormal(tempids[0], tempids[1], tempids[2]);
					assignnormal(tempids[0], tempids[2], tempids[3]);
					assignnormal(tempids[0], tempids[3], tempids[4]);
					assignnormal(tempids[0], tempids[4], tempids[5]);
					assignnormal(tempids[0], tempids[5], tempids[6]);
					assignnormal(tempids[0], tempids[6], tempids[7]);
					assignnormal(tempids[0], tempids[7], tempids[8]);
					assignnormal(tempids[0], tempids[8], tempids[1]);

					// Calculate the Sum of the normals 

					for (int i = 0; i < 8; i++) {
						fx = fx + tempnormals[i].getx();
						fy = fy + tempnormals[i].gety();
						fz = fz + tempnormals[i].getz();
					}

					fx = fx / 8;
					fy = fy / 8;
					fz = fz / 8;

					// Assign final normal to the point.
					point3d[a].normal[0] = Normal3d(fx, fy, fz);
					
				}

				else if ((h>0) && (h < gw) && (w == gw)) {
					// right column
					tempids.push_back(a);			// 0
					tempids.push_back(a + gw + 1);	// 1
					tempids.push_back(a + gw);		// 2
					tempids.push_back(a - 1);		// 3
					tempids.push_back(a - gw - 2);	// 4
					tempids.push_back(a - gw - 1);	// 5

					// Calculate normals for each of the tempids.
					// Specific assignment

					assignnormal(tempids[0], tempids[1], tempids[2]);
					assignnormal(tempids[0], tempids[2], tempids[3]);
					assignnormal(tempids[0], tempids[3], tempids[4]);
					assignnormal(tempids[0], tempids[4], tempids[5]);

					// Calculate the Sum of the normals 
					for (int i = 0; i < 4; i++) {
						fx = fx + tempnormals[i].getx();
						fy = fy + tempnormals[i].gety();
						fz = fz + tempnormals[i].getz();
					}

					fx = fx / 4;
					fy = fy / 4;
					fz = fz / 4;

					// Assign final normal to the point.
					point3d[a].normal[0] = Normal3d(fx, fy, fz);
					
				}

				else if ((h == gw) && (w == 0)) {
					// top row far left
					tempids.push_back(a);				// 0
					tempids.push_back(a - gw - 1);	// 1
					tempids.push_back(a - gw);		// 2
					tempids.push_back(a + 1);				// 3				

					// Calculate normals for each of the tempids.
					// Specific assignment

					assignnormal(tempids[0], tempids[1], tempids[2]);
					assignnormal(tempids[0], tempids[2], tempids[3]);

					// Calculate the Sum of the normals 

					for (int i = 0; i < 2; i++) {
						fx = fx + tempnormals[i].getx();
						fy = fy + tempnormals[i].gety();
						fz = fz + tempnormals[i].getz();
					}

					fx = fx / 2;
					fy = fy / 2;
					fz = fz / 2;

					// Assign final normal to the point.
					point3d[a].normal[0] = Normal3d(fx, fy, fz);
					
				}

				else if ((h == gw) && (w > 0) && (w < gw)) {
					// top row 
					tempids.push_back(a);				// 0
					tempids.push_back(a + 1);				// 1
					tempids.push_back(a - gw);		// 2
					tempids.push_back(a - gw - 1);	// 3
					tempids.push_back(a - gw - 2);	// 4
					tempids.push_back(a - 1);				// 5

					// Calculate normals for each of the tempids.
					// Specific assignment

					assignnormal(tempids[0], tempids[1], tempids[2]);
					assignnormal(tempids[0], tempids[2], tempids[3]);
					assignnormal(tempids[0], tempids[3], tempids[4]);
					assignnormal(tempids[0], tempids[4], tempids[5]);

					// Calculate the Sum of the normals 

					for (int i = 0; i < 4; i++) {
						fx = fx + tempnormals[i].getx();
						fy = fy + tempnormals[i].gety();
						fz = fz + tempnormals[i].getz();
					}

					fx = fx / 4;
					fy = fy / 4;
					fz = fz / 4;

					// Assign final normal to the point.
					point3d[a].normal[0] = Normal3d(fx, fy, fz);
					
				}

				else if ((h == gw) && (w == gw)) {
					// top row far right
					tempids.push_back(a);				// 0
					tempids.push_back(a - 1);				// 1
					tempids.push_back(a - gw - 2);	// 2
					tempids.push_back(a - gw - 1);	// 3					

					// Calculate normals for each of the tempids.
					// Specific assignment

					assignnormal(tempids[0], tempids[1], tempids[2]);
					assignnormal(tempids[0], tempids[2], tempids[3]);

					// Calculate the Sum of the normals 

					for (int i = 0; i < 2; i++) {
						fx = fx + tempnormals[i].getx();
						fy = fy + tempnormals[i].gety();
						fz = fz + tempnormals[i].getz();
					}

					fx = fx / 2;
					fy = fy / 2;
					fz = fz / 2;

					// Assign final normal to the point.
					point3d[a].normal[0] = Normal3d(fx, fy, fz);
					
				}

				tempids.clear();
				tempnormals.clear();	

				fx = 0, fy = 0; fz = 0;
				a++;
			}
		}
		// end of assignment				
	}

	Point3d normals(int v1, int v2, int v3) {		

		Point3d u = (point3d[v2] - point3d[v1]);
		Point3d v = (point3d[v3] - point3d[v1]);
				
		double nx = (u.gety() * v.getz()) - (u.getz() * v.gety());
		double ny = (u.getz() * v.getx()) - (u.getx() * v.getz());
		double nz = (u.getx() * v.gety()) - (u.gety() * v.getx());		

		double l = sqrt((nx*nx) + (ny*ny) + (nz*nz));

		nx /= l;
		ny /= l;
		nz /= l;

		return Point3d(nx, ny, nz);

	}

	void assignnormal(int v1, int v2, int v3) {		
		tempnormals.push_back(normals(v1, v2, v3));
	}	

	int draw(SpaceCamera &spacecamera) {

		int num = glGenLists(1);
		glNewList(num, GL_COMPILE);

		// the triangle selected is the referance provided by the referance ary vector.		
		glBegin(GL_TRIANGLES);
		for (int i = 0; i < triangle3d.size(); i++) {
				// v1
				//glColor3d(0.0, 0.5, 0.0);
				glNormal3d(point3d[triangle3d[i].getv1()].normal[0].getx(),
					point3d[triangle3d[i].getv1()].normal[0].gety(),
					point3d[triangle3d[i].getv1()].normal[0].getz());

				glVertex3d(point3d[triangle3d[i].getv1()].getx(),
					point3d[triangle3d[i].getv1()].gety(),
					point3d[triangle3d[i].getv1()].getz());
				glTexCoord2f(0.0, 0.0);

				// v2
				//glColor3d(0.0, 0.5, 0.0);
				glNormal3d(point3d[triangle3d[i].getv2()].normal[0].getx(),
					point3d[triangle3d[i].getv2()].normal[0].gety(),
					point3d[triangle3d[i].getv2()].normal[0].getz());

				glVertex3d(point3d[triangle3d[i].getv2()].getx(),
					point3d[triangle3d[i].getv2()].gety(),
					point3d[triangle3d[i].getv2()].getz());
				glTexCoord2f(1.0, 0.0);

				// v3
				//glColor3d(0.0, 0.5, 0.0);
				glNormal3d(point3d[triangle3d[i].getv3()].normal[0].getx(),
					point3d[triangle3d[i].getv3()].normal[0].gety(),
					point3d[triangle3d[i].getv3()].normal[0].getz());

				glVertex3d(point3d[triangle3d[i].getv3()].getx(),
					point3d[triangle3d[i].getv3()].gety(),
					point3d[triangle3d[i].getv3()].getz());
				glTexCoord2f(0.0, 1.0);
		}
		glEnd();
		glEndList();

		return num;
	}

	int drawwater(SpaceCamera &spacecamera) {		

		int num = glGenLists(1);
		glNewList(num, GL_COMPILE);

		// the triangle selected is the referance provided by the referance ary vector.		
		glBegin(GL_TRIANGLES);
		for (int i = 0; i < triangle3d.size(); i++) {
			if ((waterheight[triangle3d[i].getv1()] > 0.0) ||
				(waterheight[triangle3d[i].getv2()] > 0.0) ||
				(waterheight[triangle3d[i].getv3()] > 0.0)) {
				// v1
				//glColor3d(0.3, 0.4, 0.8);
				glVertex3d(point3d[triangle3d[i].getv1()].getx(),
					point3d[triangle3d[i].getv1()].gety() + waterheight[triangle3d[i].getv1()],
					point3d[triangle3d[i].getv1()].getz());
				glTexCoord2f(0.0, 0.0);

				// v2
				//glColor3d(0.3, 0.4, 0.8);
				glVertex3d(point3d[triangle3d[i].getv2()].getx(),
					point3d[triangle3d[i].getv2()].gety() + waterheight[triangle3d[i].getv2()],
					point3d[triangle3d[i].getv2()].getz());
				glTexCoord2f(1.0, 0.0);

				// v3
				//glColor3d(0.3, 0.4, 0.8);
				glVertex3d(point3d[triangle3d[i].getv3()].getx(),
					point3d[triangle3d[i].getv3()].gety() + waterheight[triangle3d[i].getv3()],
					point3d[triangle3d[i].getv3()].getz());
				glTexCoord2f(0.0, 1.0);
			}
		}
		glEnd();
		glEndList();

		return num;
	}

};

#endif