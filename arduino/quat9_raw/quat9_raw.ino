#include "ICM_20948.h" // Click here to get the library: http://librarymanager/All#SparkFun_ICM_20948_IMU

//#define AD0_VAL 1      // The value of the last bit of the I2C address.

ICM_20948_I2C myICM; // Otherwise create an ICM_20948_I2C object

// gravity sucks
//long acl_offset_x = 0;
//long acl_offset_y = 0;
//long acl_offset_z = 0;
long gyr_th_x = 0;
long gyr_th_y = 0;
long gyr_th_z = 0;
int counter = 0;
void setup() {

  Serial.begin(38400); // Start the serial console // 115200
  while (!Serial);

  //delay(100);

  Wire.begin();
  Wire.setClock(400000);

  //bool initialized = false;
  //while (!initialized) {
  //while (1) {
  //  // Initialize the ICM-20948
  //  // If the DMP is enabled, .begin performs a minimal startup. We need to configure the sample mode etc. manually.
  //  myICM.begin(Wire, 1); // AD0_VAL 

  //  //Serial.print("Init: ");
  //  //Serial.println(myICM.statusString());

  //  //if (myICM.status != ICM_20948_Stat_Ok) {
  //  //  //Serial.println("init fail");
  //  //  delay(500);
  //  //} else {
  //  //  break;
  //  //  //initialized = true;
  //  //}

  //  if (myICM.status == ICM_20948_Stat_Ok) {
  //     break;
  //  }
  //  //delay(500);
  //} // !initialized

  // initialize
  do {
     myICM.begin(Wire, 1); // AD0_VAL = 1
  } while (myICM.status != ICM_20948_Stat_Ok);

  //Serial.println(F("Device connected!"));

  //bool success = true; // Use success to show if the DMP configuration was successful

  // Initialize the DMP. initializeDMP is a weak function. You can overwrite it if you want to e.g. to change the sample rate
  bool success = (myICM.initializeDMP() == ICM_20948_Stat_Ok);
  //myICM.initializeDMP(); // minimal
  // DMP sensor options are defined in ICM_20948_DMP.h
  //    INV_ICM20948_SENSOR_ACCELEROMETER               (16-bit accel)
  //    INV_ICM20948_SENSOR_GYROSCOPE                   (16-bit gyro + 32-bit calibrated gyro)
  //    INV_ICM20948_SENSOR_RAW_ACCELEROMETER           (16-bit accel)
  //    INV_ICM20948_SENSOR_RAW_GYROSCOPE               (16-bit gyro + 32-bit calibrated gyro)
  //    INV_ICM20948_SENSOR_MAGNETIC_FIELD_UNCALIBRATED (16-bit compass)
  //    INV_ICM20948_SENSOR_GYROSCOPE_UNCALIBRATED      (16-bit gyro)
  //    INV_ICM20948_SENSOR_STEP_DETECTOR               (Pedometer Step Detector)
  //    INV_ICM20948_SENSOR_STEP_COUNTER                (Pedometer Step Detector)
  //    INV_ICM20948_SENSOR_GAME_ROTATION_VECTOR        (32-bit 6-axis quaternion)
  //    INV_ICM20948_SENSOR_ROTATION_VECTOR             (32-bit 9-axis quaternion + heading accuracy)
  //    INV_ICM20948_SENSOR_GEOMAGNETIC_ROTATION_VECTOR (32-bit Geomag RV + heading accuracy)
  //    INV_ICM20948_SENSOR_GEOMAGNETIC_FIELD           (32-bit calibrated compass)
  //    INV_ICM20948_SENSOR_GRAVITY                     (32-bit 6-axis quaternion)
  //    INV_ICM20948_SENSOR_LINEAR_ACCELERATION         (16-bit accel + 32-bit 6-axis quaternion)
  //    INV_ICM20948_SENSOR_ORIENTATION                 (32-bit 9-axis quaternion + heading accuracy)

  // Enable the DMP orientation sensor
  success &= (myICM.enableDMPSensor(INV_ICM20948_SENSOR_ORIENTATION) == ICM_20948_Stat_Ok);
  //myICM.enableDMPSensor(INV_ICM20948_SENSOR_ORIENTATION); minimal
  // Enable any additional sensors / features
  //success &= (myICM.enableDMPSensor(INV_ICM20948_SENSOR_RAW_GYROSCOPE) == ICM_20948_Stat_Ok);
  //success &= (myICM.enableDMPSensor(INV_ICM20948_SENSOR_RAW_ACCELEROMETER) == ICM_20948_Stat_Ok);
  //success &= (myICM.enableDMPSensor(INV_ICM20948_SENSOR_MAGNETIC_FIELD_UNCALIBRATED) == ICM_20948_Stat_Ok);

  // Configuring DMP to output data at multiple ODRs:
  // DMP is capable of outputting multiple sensor data at different rates to FIFO.
  // Setting value can be calculated as follows:
  // Value = (DMP running rate / ODR ) - 1
  // E.g. For a 5Hz ODR rate when DMP is running at 55Hz, value = (55/5) - 1 = 10.
  success &= (myICM.setDMPODRrate(DMP_ODR_Reg_Quat9, 0) == ICM_20948_Stat_Ok); // Set to the maximum
  //myICM.setDMPODRrate(DMP_ODR_Reg_Quat9, 0); minimal
  //success &= (myICM.setDMPODRrate(DMP_ODR_Reg_Accel, 0) == ICM_20948_Stat_Ok); // Set to the maximum
  //success &= (myICM.setDMPODRrate(DMP_ODR_Reg_Gyro, 0) == ICM_20948_Stat_Ok); // Set to the maximum
  //success &= (myICM.setDMPODRrate(DMP_ODR_Reg_Gyro_Calibr, 0) == ICM_20948_Stat_Ok); // Set to the maximum
  //success &= (myICM.setDMPODRrate(DMP_ODR_Reg_Cpass, 0) == ICM_20948_Stat_Ok); // Set to the maximum
  //success &= (myICM.setDMPODRrate(DMP_ODR_Reg_Cpass_Calibr, 0) == ICM_20948_Stat_Ok); // Set to the maximum

  // Enable the FIFO
  success &= (myICM.enableFIFO() == ICM_20948_Stat_Ok);
  //myICM.enableFIFO(); // minimal
  // Enable the DMP
  success &= (myICM.enableDMP() == ICM_20948_Stat_Ok);
  //myICM.enableDMP(); // minimal
  // Reset DMP
  success &= (myICM.resetDMP() == ICM_20948_Stat_Ok);
  //myICM.resetDMP(); // minimal
  // Reset FIFO
  success &= (myICM.resetFIFO() == ICM_20948_Stat_Ok);
  //myICM.resetFIFO(); // minimal
  // Check success
  //if (success) {
  //  Serial.println("DMP en");
  //} else {
  if (!success) {
     Serial.println("DMP fail");
  }
  while(!success); // failed initialization
  //if (!success) {
  //  //Serial.println("DMP fail");
  //  //Serial.println(F("Please check that you have uncommented line 29 (#define ICM_20948_USE_DMP) in ICM_20948_C.h..."));
  //  while (1); // Do nothing more
  //}

  // gravity sucks
  //for(int i = 0; i < 10; i++) {
  //   while (!myICM.dataReady());
  //   myICM.getAGMT(); // poll IMU
  //   acl_offset_x += myICM.agmt.acc.axes.x;
  //   acl_offset_y += myICM.agmt.acc.axes.y;
  //   acl_offset_z += myICM.agmt.acc.axes.z;
  //}
  //acl_offset_x /= 10;
  //acl_offset_y /= 10;
  //acl_offset_z /= 10;
  
  //long gyr_offset_x = 0;
  //long gyr_offset_y = 0; 
  //for(int i = 0; i < 10; i++) {
  //   while (!myICM.dataReady());
  //   myICM.getAGMT(); // poll IMU
  //   gyr_offset_x += myICM.agmt.gyr.axes.x;
  //   gyr_offset_y += myICM.agmt.gyr.axes.y;
  //   //gyr_offset_z += myICM.agmt.gyr.axes.z;
  //}
  //gyr_offset_x /= 10;
  //gyr_offset_y /= 10;
  //gyr_th_x = abs(gyr_offset_x) + 1000;
  //gyr_th_y = abs(gyr_offset_y) + 1000; 
  //gyr_th = sqrt(gyr_offset_x/131 * gyr_offset_x/131 + gyr_offset_y/131 * gyr_offset_y/131);

  gyr_th_x = 1000;
  gyr_th_y = 1000;
  gyr_th_z = 1000;

  Serial.print("gyr threshold: ");
  Serial.print(gyr_th_x);
  Serial.print(" , ");
  Serial.print(gyr_th_y);
  Serial.print(" , ");
  Serial.println(gyr_th_z);
} // endfunction setup

void loop() {
   icm_20948_DMP_data_t data;
   myICM.readDMPdataFromFIFO(&data);

   if ((myICM.status == ICM_20948_Stat_Ok) || (myICM.status == ICM_20948_Stat_FIFOMoreDataAvail)) { // Was valid data available?

      if ((data.header & DMP_header_bitmap_Quat9) > 0) { 

         if (myICM.dataReady()) {
            myICM.getAGMT(); // poll IMU

            if ((abs(myICM.agmt.gyr.axes.x) > gyr_th_x) || (abs(myICM.agmt.gyr.axes.y) > gyr_th_y) || (abs(myICM.agmt.gyr.axes.z) > gyr_th_z)) {
               if (counter > 10) {
                  Serial.print("#!"); Serial.print(",");
                  Serial.print("zero-counter: "); Serial.println(counter);
               }
               Serial.print(data.Quat9.Data.Q1); Serial.print(",");
               Serial.print(data.Quat9.Data.Q2); Serial.print(",");
               Serial.println(data.Quat9.Data.Q3); // dont use println if want to print raw values too
               counter = 0;
            } else {
               counter++;
            }
         } // myICM.dataReady
      } // (data.header & DMP_header_bitmap_Quat9) > 0
   } // Was valid data available?
} // loop()
