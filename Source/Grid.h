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

	double gradient;
	
public:
	Point3d(double nx, double ny, double nz) : x(nx), y(ny), z(nz) {
		gradient=0;
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

	Point3d operator/(double a) {
		Point3d b;

		b.x = this->x / a;
		b.y = this->y / a;
		b.z = this->z / a;
		return b;
	}

	double getgradient() { return gradient; }	
	void setgradient(double ng) { gradient = ng; }

	Normal3d normal[1];

	struct {
		double r;
		double g;
		double b;

		void setR(double nr) {
			r=nr;
		}
		void setG(double ng) {
			g=ng;
		}
		void setB(double nb) {
			b=nb;
		}

		double getR() { return r; }
		double getG() { return g; }
		double getB() { return b; }
	} GradientColour;
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

		switch (rangevalues) {

			case 0:
				for (int i = 1; i < length; i++) {
					double newvalue_x = xvalues[i - 1] + rn.GenerateUniformRealDistribution(0.0, 0.025);					
					xvalues.push_back(newvalue_x);

					double newvalue_y = yvalues[i - 1] + rn.GenerateUniformRealDistribution(0.0, 0.025);			
					yvalues.push_back(newvalue_y);
				}
				break;

			case 1:
				for (int i = 1; i < length; i++) {
					
					double newvalue_x = xvalues[i - 1] + rn.GenerateUniformRealDistribution(-0.025, 0.0);
					xvalues.push_back(newvalue_x);
					
					double newvalue_y = yvalues[i - 1] + rn.GenerateUniformRealDistribution(0.0, 0.025);
					yvalues.push_back(newvalue_y);
				}
				break;

			case 2:
				for (int i = 1; i < length; i++) {
					
					double newvalue_x = xvalues[i - 1] + rn.GenerateUniformRealDistribution(-0.025, 0.0);
					xvalues.push_back(newvalue_x);
										
					double newvalue_y = yvalues[i - 1] + rn.GenerateUniformRealDistribution(-0.025, 0.0);
					yvalues.push_back(newvalue_y);
				}
				break;

			case 3:
				for (int i = 1; i < length; i++) {
					
					double newvalue_x = xvalues[i - 1] + rn.GenerateUniformRealDistribution(0.0, 0.025);
					xvalues.push_back(newvalue_x);
										
					double newvalue_y = yvalues[i - 1] + rn.GenerateUniformRealDistribution(-0.025, 0.0);
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

		double offset = rn.GenerateUniformIntDistribution(0, 12);

		for (int i = 0; i < zvalues.size(); i++) {

			double noise = amplitude * sin((getZValue(i) * frequancy) + offset);

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

	RandomNumbers rn;
	double frequancy;
	double amplitude;
		

public:
	CGrid(int nlength) : length(nlength) {

		frequancy = 1.0;
		amplitude = 1.0;		

		initializegrid();
		
	}

	~CGrid() {}

	std::vector<Point3d> point3d;
	std::vector<Triangle3d> triangle3d;	

	std::vector<NoiseContainer> noisecontainer;
	std::vector<double> sumofnoise;

	std::vector<int> tempids;
	std::vector<Point3d> tempnormals;
	Point3d resultantnormal;
	std::vector<double> tempgradients;

	std::vector<double> noise;

	
	void addPoint3d(double nx, double ny, double nz) { point3d.push_back(Point3d(nx, ny, nz)); }
	void addTriangle3d(int nv1, int nv2, int nv3)    { triangle3d.push_back(Triangle3d(nv1, nv2, nv3)); }

	void initializegrid() {
		
		for (int i = 0; i < length * length; i++) {

			sumofnoise.push_back(0.0);

		}

		double max = 0;
		double min = 0;

		for (int i = 0; i < 8; i++) {
			
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


		/*for (int i = 0; i < sumofnoise[i]; i++) {

			if (sumofnoise[i] > max) { max = sumofnoise[i]; }
			if (sumofnoise[i] < min) { min = sumofnoise[i]; }

		}
		for (int i = 0; i < sumofnoise[i]; i++) {

			sumofnoise[i] = (sumofnoise[i] - min) / (max - min);

		}*/

		for (int h = 0; h < length; h++) {

			for (int w = 0; w < length; w++) {

				addPoint3d(w, 0.0, h);

			}

		}				
		
		
		for (int p = 0; p < 8; p++) {
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
		
		for (int i = 0; i < sumofnoise[i]; i++) {

			if (sumofnoise[i] > max) { max = sumofnoise[i]; }
			if (sumofnoise[i] < min) { min = sumofnoise[i]; }

		}
		for (int i = 0; i < sumofnoise[i]; i++) {

			sumofnoise[i] = (sumofnoise[i] - min) / (max - min);

		}

		for (int i = 0; i < length * length; i++) {			

			point3d[i] = Point3d(point3d[i].getx(), sumofnoise[i], point3d[i].getz());

			point3d[i] = Point3d(point3d[i].getx(), 50*pow(point3d[i].gety(), 2.0), point3d[i].getz());

			//point3d[i] = Point3d(point3d[i].getx(), 10.0 * sumofnoise[i], point3d[i].getz());			

		}

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

	void vertexnormals() {
		int gw = length - 1;
		int a = 0;
		for (int h = 0; h < length; h++) {
			for (int w = 0; w < length; w++) {

				//std::vector<int> tempids;
				//std::vector<Point3d> tempnormals;

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

					calculategradient(1);
					calculategradient(2);
					calculategradient(3);

					double totalgradient = 0.0;

					for (int i = 0; i < tempgradients.size(); i++) {
						totalgradient = totalgradient + tempgradients[i];
					}

					point3d[a].setgradient(totalgradient);
					point3d[a].GradientColour.setR(totalgradient);
					point3d[a].GradientColour.setG(0.5);
					point3d[a].GradientColour.setB(0.5);

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

					calculategradient(1);
					calculategradient(2);
					calculategradient(3);
					calculategradient(4);
					calculategradient(5);

					double totalgradient = 0.0;

					for (int i = 0; i < tempgradients.size(); i++) {
						totalgradient = totalgradient + tempgradients[i];
					}

					point3d[a].setgradient(totalgradient);

					point3d[a].GradientColour.setR(totalgradient);
					point3d[a].GradientColour.setG(0.5);
					point3d[a].GradientColour.setB(0.5);
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

					calculategradient(1);
					calculategradient(2);
					calculategradient(3);

					double totalgradient = 0.0;

					for (int i = 0; i < tempgradients.size(); i++) {
						totalgradient = totalgradient + tempgradients[i];
					}

					point3d[a].setgradient(totalgradient);

					point3d[a].GradientColour.setR(totalgradient);
					point3d[a].GradientColour.setG(0.5);
					point3d[a].GradientColour.setB(0.5);
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

					calculategradient(1);
					calculategradient(2);
					calculategradient(3);
					calculategradient(4);
					calculategradient(5);

					double totalgradient = 0.0;

					for (int i = 0; i < tempgradients.size(); i++) {
						totalgradient = totalgradient + tempgradients[i];
					}

					point3d[a].setgradient(totalgradient);

					point3d[a].GradientColour.setR(totalgradient);
					point3d[a].GradientColour.setG(0.5);
					point3d[a].GradientColour.setB(0.5);


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

					calculategradient(1);
					calculategradient(2);
					calculategradient(3);
					calculategradient(4);
					calculategradient(5);
					calculategradient(6);
					calculategradient(7);
					calculategradient(8);

					double totalgradient = 0.0;

					for (int i = 0; i < tempgradients.size(); i++) {
						totalgradient = totalgradient + tempgradients[i];
					}

					point3d[a].setgradient(totalgradient);

					point3d[a].GradientColour.setR(totalgradient);
					point3d[a].GradientColour.setG(0.5);
					point3d[a].GradientColour.setB(0.5);
				}
				else if ((h>0) && (h < gw) && (w == gw)) {
					// right column
					tempids.push_back(a);				// 0
					tempids.push_back(a + gw + 1);	// 1
					tempids.push_back(a + gw);		// 2
					tempids.push_back(a - 1);				// 3
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

					calculategradient(1);
					calculategradient(2);
					calculategradient(3);
					calculategradient(4);
					calculategradient(5);

					double totalgradient = 0.0;

					for (int i = 0; i < tempgradients.size(); i++) {
						totalgradient = totalgradient + tempgradients[i];
					}

					point3d[a].setgradient(totalgradient);

					point3d[a].GradientColour.setR(totalgradient);
					point3d[a].GradientColour.setG(0.5);
					point3d[a].GradientColour.setB(0.5);
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

					calculategradient(1);
					calculategradient(2);
					calculategradient(3);

					double totalgradient = 0.0;

					for (int i = 0; i < tempgradients.size(); i++) {
						totalgradient = totalgradient + tempgradients[i];
					}

					point3d[a].setgradient(totalgradient);

					point3d[a].GradientColour.setR(totalgradient);
					point3d[a].GradientColour.setG(0.5);
					point3d[a].GradientColour.setB(0.5);
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

					calculategradient(1);
					calculategradient(2);
					calculategradient(3);
					calculategradient(4);
					calculategradient(5);

					double totalgradient = 0.0;

					for (int i = 0; i < tempgradients.size(); i++) {
						totalgradient = totalgradient + tempgradients[i];
					}

					point3d[a].setgradient(totalgradient);

					point3d[a].GradientColour.setR(totalgradient);
					point3d[a].GradientColour.setG(0.5);
					point3d[a].GradientColour.setB(0.5);
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

					calculategradient(1);
					calculategradient(2);
					calculategradient(3);

					double totalgradient = 0.0;

					for (int i = 0; i < tempgradients.size(); i++) {
						totalgradient = totalgradient + tempgradients[i];
					}

					point3d[a].setgradient(totalgradient);

					point3d[a].GradientColour.setR(totalgradient);
					point3d[a].GradientColour.setG(0.5);
					point3d[a].GradientColour.setB(0.5);
				}

				tempids.clear();
				tempnormals.clear();
				tempgradients.clear();

				fx = 0, fy = 0; fz = 0;
				a++;
			}
		}
		// end of assignment				
	}

	void normals(int v1, int v2, int v3) {

		// v1 centre point
		// v2 lhs
		// v3 rhs

		// The maths here needs checking. Although this method is correct 
		// it might produce unintended results if data is passed incorretly

		// example where v1 10, v2 20, v3 20
		// u = v2-v1  
		// u = 10

		// v = v3-v1
		// v = 10

		// if v1 are 20 and v2 10
		// u = -10
		// v = -10

		// dispite being the same vector?

		// u = (v2+v1) /2 
		// v = (v3+v1) /2 

		// solves this?

		double u_x = point3d[v2].getx() - point3d[v1].getx();
		double u_y = point3d[v2].gety() - point3d[v1].gety();
		double u_z = point3d[v2].getz() - point3d[v1].getz();

		double v_x = point3d[v3].getx() - point3d[v1].getx();
		double v_y = point3d[v3].gety() - point3d[v1].gety();
		double v_z = point3d[v3].getz() - point3d[v1].getz();

		double nx = (u_y*v_z) - (u_z*v_y);
		double ny = (u_z*v_x) - (u_x*v_z);
		double nz = (u_x*v_y) - (u_y*v_x);

		double l = sqrt((nx*nx) + (ny*ny) + (nz*nz));

		nx /= l;
		ny /= l;
		nz /= l;

		resultantnormal = Point3d(nx, ny, nz);

	}

	void assignnormal(int v1, int v2, int v3) {
		normals(v1, v2, v3);
		tempnormals.push_back(resultantnormal);
	}

	void calculategradient(int v1) {
		double dy = point3d[0].gety() - point3d[v1].gety();
		tempgradients.push_back(dy);
	}

	int draw(SpaceCamera &spacecamera) {

		int num = glGenLists(1);
		glNewList(num, GL_COMPILE);

		// the triangle selected is the referance provided by the referance ary vector.		
		glBegin(GL_TRIANGLES);
		for (int i = 0; i < triangle3d.size(); i++) {

			// v1
			glNormal3d(point3d[triangle3d[i].getv1()].normal[0].getx(),
				point3d[triangle3d[i].getv1()].normal[0].gety(),
				point3d[triangle3d[i].getv1()].normal[0].getz());

			glVertex3d(point3d[triangle3d[i].getv1()].getx(),
				point3d[triangle3d[i].getv1()].gety(),
				point3d[triangle3d[i].getv1()].getz());
			glTexCoord2f(0.0, 0.0);

			// v2
			glNormal3d(point3d[triangle3d[i].getv2()].normal[0].getx(),
				point3d[triangle3d[i].getv2()].normal[0].gety(),
				point3d[triangle3d[i].getv2()].normal[0].getz());

			glVertex3d(point3d[triangle3d[i].getv2()].getx(),
				point3d[triangle3d[i].getv2()].gety(),
				point3d[triangle3d[i].getv2()].getz());
			glTexCoord2f(1.0, 0.0);

			// v3
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

};

#endif
