#include "ICM_20948.h" // Click here to get the library: http://librarymanager/All#SparkFun_ICM_20948_IMU

//#define AD0_VAL 1      // The value of the last bit of the I2C address.

ICM_20948_I2C myICM; // Otherwise create an ICM_20948_I2C object

long gyr_th_x = 0;
long gyr_th_y = 0;
long gyr_th_z = 0;
int counter = 0;
void setup() {

  Serial.begin(38400); 
  while (!Serial);

  Wire.begin();
  Wire.setClock(400000);

  // initialize
  do {
     myICM.begin(Wire, 1); // AD0_VAL = 1
  } while (myICM.status != ICM_20948_Stat_Ok);

  // Initialize the DMP. initializeDMP is a weak function. You can overwrite it if you want to e.g. to change the sample rate
  bool success = (myICM.initializeDMP() == ICM_20948_Stat_Ok);
  //    INV_ICM20948_SENSOR_ORIENTATION                 (32-bit 9-axis quaternion + heading accuracy)

  // Enable the DMP orientation sensor
  success &= (myICM.enableDMPSensor(INV_ICM20948_SENSOR_ORIENTATION) == ICM_20948_Stat_Ok);
  // Configuring DMP to output data at multiple ODRs:
  // DMP is capable of outputting multiple sensor data at different rates to FIFO.
  // Setting value can be calculated as follows:
  // Value = (DMP running rate / ODR ) - 1
  // E.g. For a 5Hz ODR rate when DMP is running at 55Hz, value = (55/5) - 1 = 10.
  success &= (myICM.setDMPODRrate(DMP_ODR_Reg_Quat9, 0) == ICM_20948_Stat_Ok); // Set to maximum ODR
  // Enable the FIFO
  success &= (myICM.enableFIFO() == ICM_20948_Stat_Ok);
  // Enable the DMP
  success &= (myICM.enableDMP() == ICM_20948_Stat_Ok);
  // Reset DMP
  success &= (myICM.resetDMP() == ICM_20948_Stat_Ok);
  // Reset FIFO
  success &= (myICM.resetFIFO() == ICM_20948_Stat_Ok);
  if (!success) {
     Serial.println("DMP fail");
  }
  while(!success); // failed initialization
  
  gyr_th_x = 1000;
  gyr_th_y = 1000;
  gyr_th_z = 1000;

  Serial.println("x,y,z");

  //Serial.print("gyr threshold: ");
  //Serial.print(gyr_th_x);
  //Serial.print(" , ");
  //Serial.print(gyr_th_y);
  //Serial.print(" , ");
  //Serial.println(gyr_th_z);
} // endfunction setup

void loop() {
   icm_20948_DMP_data_t data;
   myICM.readDMPdataFromFIFO(&data);

   if ((myICM.status == ICM_20948_Stat_Ok) || (myICM.status == ICM_20948_Stat_FIFOMoreDataAvail)) { // Was valid data available?

      if ((data.header & DMP_header_bitmap_Quat9) > 0) { 

         if (myICM.dataReady()) {
            myICM.getAGMT(); // poll IMU

            if ((abs(myICM.agmt.gyr.axes.x) > gyr_th_x) || (abs(myICM.agmt.gyr.axes.y) > gyr_th_y) || (abs(myICM.agmt.gyr.axes.z) > gyr_th_z)) {
               //Serial.print(myICM.agmt.acc.axes.x); Serial.print(",");
              // Serial.print(myICM.agmt.acc.axes.y); Serial.print(",");
               //Serial.println(myICM.agmt.acc.axes.z);
               if (counter > 10) {
                  Serial.print("#!"); Serial.print(",");
                  Serial.print("zero-counter: "); Serial.println(counter);
               }
               //if (abs(myICM.agmt.gyr.axes.z) > 8000) {
               //   Serial.print("#!,gyr_z: ");
               //   Serial.println(myICM.agmt.gyr.axes.z);
               //}
               //if (myICM.agmt.gyr.axes.z > 6000) {
                  //Serial.println("0,0,0");
               //} //else {
                  Serial.print(data.Quat9.Data.Q1); Serial.print(",");
                  Serial.print(data.Quat9.Data.Q2); Serial.print(",");
                  Serial.print(data.Quat9.Data.Q3); // dont use println if want to print raw values too
                  Serial.print(","); Serial.println(myICM.agmt.acc.axes.x);
               //}
               counter = 0;
            } else {
               counter++;
            }
            //
         } // myICM.dataReady
      } // (data.header & DMP_header_bitmap_Quat9) > 0
   } // Was valid data available?
} // loop()
