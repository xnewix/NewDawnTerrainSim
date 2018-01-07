#include <iostream>
#include <iomanip>
#include <vector>
#include <random>

const double PI = 3.1415926535897;

class Vec3 {
private:
  double x, y, z;
public:
  Vec3(double nx, double ny, double nz) : x(nx), y(ny), z(nz) {    
  };
  Vec3();
  //~Vec3();
  
  double getX() { return x; }
  double getY() { return y; }
  double getZ() { return z; }

};

class CGrid {
private:
  int length;
  double latitude;
public:
    CGrid(int nlength) : length(nlength) {    
    };
    CGrid();
    //~CGrid();
    std::vector<Vec3>   points;
    std::vector<double> temperature;
    std::vector<double> waterheight;
    std::vector<double> sin_percentage;
    std::vector<double> moisture;
    
    void init() {
        latitude = 45.0;
        double latitudebase_temp = latitude * 0.5;
      
        int m = 0;
        for(int h = 0; h < length; h++) {
            for(int w = 0; w < length; w++) {
                points.push_back(Vec3(w, w * h, h));
                temperature.push_back((latitudebase_temp - (h * 0.0000045045)) - 
                                      (points[points.size()-1].getY() * 0.0098));
                 
                // In the Northen hemisphere the latitude moves futher from 0 degrees, 
                // the temperature should decrease
              
                // 0.5 degrees increase per degree of latitude?
                // 111 km in a point of latitude.  
                // therefore 0.0000045045 degrees temperature change per meter.
                
                // 45 lat 25 degrees C 
              
                //the X axis auto
                /*
                std::cout << std::setprecision(16)
                          << "  x: " << points[points.size()-1].getX() 
                          << "  y: " << points[points.size()-1].getY()
                          << "  z: " << points[points.size()-1].getZ() 
                          << "  c: " << temperature[points.size()-1] 
                          << std::endl;
                */
                moisture.push_back(0);
                waterheight.push_back(5);
                m++;
            }
        }
        
        // Simulate Yearly Temperature Cycle
        for(double i = 0.0; i < 365.0; i++) {
            // generate the percentage table 
          
            double fraction = i / 365.0;
            double percentage = fraction * (2.0 * PI);
            double sinvalue = sin(percentage);
            double sinlowerrange;
            //std::cout << "i: " << i << ", " << fraction << ", " << percentage << ", " << sinvalue << std::endl;
          
            // Calculate a flat temperature change to apply to whole map.
            // Regardless of local temerature values. 
            // Need the average centre value. Not the final temperature 
            // but the value that is applied to the base temeprature 
            
            double lat = 55.0;
            double lat_inverse = 90.0 - lat;
            double testtemperature = ((lat_inverse * 0.5) - (0 * 0.0000045045)) - (0.0 * 0.0098);
            double dailytemperature = testtemperature + (testtemperature * sinvalue);
            //std::cout << "it: " << i << ", " << testtemperature << ", " << dailytemperature <<  std::endl;
            //if(sinvalue < 0.0) { sinlowerrange = sinvalue;// + (sinvalue * 0.5); 
            //                     std::cout << "lower: " << sinlowerrange << std::endl; }
          
            // We can adjust the range of temperature values cycled 
            if(i < (365.0 / 2.0)) {
                sinlowerrange = testtemperature + (testtemperature * (sinvalue - (sinvalue * 0.5)));
                //std::cout << i << "lower: " << sinlowerrange << std::endl;
            }          
            sin_percentage.push_back(sinvalue);
        }
      
        
    }
    
    void updateTemperature(int dayofyear) {
      
        double lat = 55.0;
        double lat_inverse = 90.0 - lat;
        
       
        for(int i = 0; i < points.size(); i++) {
            double base_temperature = ((lat_inverse * 0.5) - (points[i].getZ() * 0.0000045045)) - (0.0 * 0.0098);
            double regional_temperature = lat_inverse * 0.5;
            //double dailytemperature = testtemperature + (testtemperature * sinvalue);
      
            //std::cout << temperature.size() << std::endl;
                
            temperature[i] =  base_temperature + (regional_temperature * sin_percentage[dayofyear]); 
            //std::cout << "Temperature:          " << temperature[i] << std::endl;               
            //std::cout << "Regional Temperature: " << regional_temperature << std::endl; 
            //std::cout << "Base Temperature:     " << base_temperature << std::endl; 
        }
        
    }
   
    void moveMoisture() {
        // Add Moisture to the air
        // Remove Moisture from the air
        // Move Moisture 
        std::vector<double> tmp_moisture;
        for(int h = 0; h < length; h++) {            
              
            tmp_moisture.push_back(moisture[getIndex(h, length -1)]);             
            
          
            for(int w = 0; w < length -1; w++) {
                tmp_moisture.push_back(moisture[getIndex(h, w)]);                 
            }
            // Move moisture 
            //tmp_moisture.clear();
        }
      
        for(int i = 0; i < length * length; i++) {
            
            //std::cout << tmp_moisture[i] << " : " << moisture[i] << std::endl;
        }
        moisture = tmp_moisture;
        
    }
    
    void evaporation() {
        
        for(int i = 0; i < length * length; i++) {
        
            // can only evaporate if RH < 100 % 
            double SVP = 6.11 * (pow(10.0, (7.5 * temperature[i]) / (temperature[i] + 237.3)));
        
            // or if there is water present
      
            double Tm = temperature[i] + (0.006 * points[i].getY());
            double TminusTd = (0.0023 * points[i].getY()) + (0.37 * points[i].getY()) + (0.53 * points[i].getY()) - 10.9;    
            double Alatitude = (90.0 - 25.0) + (points[i].getZ() * 0.0000045045);
            double vapour = ((((700.0 * Tm) / (100 - Alatitude)) + (15 * TminusTd)) / (80.0 - temperature[i]) / 36.0);
            double Pmbar  = 0.09806648572 * vapour;
          
            if(((vapour / 1000.0) <= waterheight[i]) && (moisture[i] + Pmbar <= SVP)) {
                moisture[i] = moisture[i] + Pmbar;
                waterheight[i] = waterheight[i] - (vapour / 1000.0);
            }
            //std::cout << temperature[i] << " : " << vapour << " : " << Pmbar << " : " << SVP << " : " << waterheight[i] << " : " << moisture[i] << std::endl;
      
            // vapour is in mm need about 220 mm to 20 pmbar = 0.2 meters 6mm a day 
        }
        std::cout << temperature[0] << " : " << waterheight[0] << " : " << moisture[0] << std::endl;
      
    }
  
   
    int getIndex(int h, int w) {
      
        int index = (h * length) + w;
        return index;
    }
    
    // represents 1 dau of the weather system 
    void cycle(int d) {
      updateTemperature(d);
      for(int i = 0; i < 360.0; i ++) {
          moveMoisture();
          evaporation();
      }
    }
    
    void ycycle() {
        for(int i = 0; i < 365.0; i++) {
            cycle(i);
        }
    }
   // Temperature should increase / decrease slighty in a cycle like motion. 
   // A seperate Air - Temperature? 
   
   // Percentage change / averaged
   // maximum of 10% averaged 6 + 8 / 2 = 7 
   // alter base value by variance +/- 
   // within limits
   
  
  
  
  
   // Cooling air is the basic method for removing moisture from the air. 
   // Warm air can hold more moistur ethan cool air. 
   // Moisture is removed from the air when it cools below its dew point
   
   // Saturated vapor pressure (SVP) = 6.11 x 10 ^ ((7.5T)/(T+237.3))
  
   // E0 = 700Tm / (100 - A) + 15(T-Td) over (80 - T)
   // Tm = T + 0.006h where h = elevation in meters
   // T-Td = 0.0023h + 0.37T + 0.53T + 0.35T - 10.9
   // A = latitude in degrees
  
   // 0.09806648572 × PmmH2O = Pmbar
  
};


CGrid cgrid(6);


int main() {
  
  
  cgrid.init();
  //cgrid.updateTemperature(0);
  //cgrid.moveMoisture();
  //cgrid.moveMoisture();
  //cgrid.evaporation();
  //cgrid.calculatetemperature(0);
  
  cgrid.cycle(0);
  std::cout << "safsffaf";
  
  
  
  
  return 0;
}
