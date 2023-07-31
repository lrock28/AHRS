# Attitude and Heading Reference System (AHRS)

   * [License](LICENSE)
   * [Changelog](CHANGELOG)

# Description

# Installation

## Matlab
Simply clone or download and the library into your  MATLAB  path folder. MATLAB 2023a was used. Include the dataset if desired for sample data.

# Constants
The following constants are defined from WGS84:

| Description | Variable | Value |
| --- | --- | --- |
| Semi-major axis, m | double SEMI_MAJOR_AXIS_LENGTH_M | 6378137.0 |
| Semi-minor axis, m | double SEMI_MINOR_AXIS_LENGTH_M | 6356752.3142 |
| Flattening | double FLATTENING | 1.0 / 298.257223563 |
| First eccentricity | double ECC | 8.1819190842622e-2 |
| First eccentricity squared | double ECC2 | 6.69437999014e-3 |
| Angular velocity of the Earth, rad/s | WE_RADPS | 7292115.0e-11 |
| Earth's gravitational constant, m^3/s^2 | GM_M3PS2 | 3986004.418e8 |

The following constants were also used:

| Description | Variable | Value |
| --- | --- | --- |
| Time Step, seconds | dt_s | 0.01 |
| Magnetic Declination, deg | magvar_deg | -3.26581 |

# Filters

## 6 State EKF AHRS

# Transformations


## Attitude
These functions convert between:
   * Euler angles ('ZYX' order rotation angles, i.e. yaw, pitch, roll)
   * Direction cosine matrix
   * Quaternions (i, j, k, w)
   * Magnetometer data
   * Accelerometer data

**bfs_euler2quaternion(const T rot1, const T rot2, const T rot3)** Converts a given set of Euler angles to quaternion. Input is the Euler angle vector (order 'ZYX'). Angle input is in radians.

```MATLAB
yaw = 0.7854; 
pitch = 0.1; 
roll = 0;
eul = [yaw, pitch, roll];
quaternion = bfs_euler2quaternion (eul);
```

**bfs_quat2eul(quaternion)** Converts a quaternion to Euler angles. Input is the quaternion. Angles are output in radians.

```MATLAB
q.x = 0.215509;
q.y = 0.432574;
q.z = 0.0846792;
q.w = 0.871358;
quaternion = [q.x, q.y, q.z, q.w];
euler = bfs_quaternion2euler(quaternion);
```

**bfs_quat2dcm(quaternion)** Convert quaternion to direction cosine matrix. Input is the quaternion and output is the direction cosine matrix.

```MATLAB
q.x = 0.215509;
q.y = 0.432574;
q.z = 0.0846792;
q.w = 0.871358;
quaternion = [q.x, q.y, q.z, q.w];
dcm = bfs_quaternion2dcm(quaternion);
```

**tilt(imu)** Convert accelerometer measurements into roll and pitch. Input is accelerometer reading in meters per second. Output is in radians.

```MATLAB
accelmps = [0.1, 0.1, 0.1];
[roll, pitch] = tilt(accelmps);
```

**psi_mag = tilt_comp (hb, phi, theta, magvar_rad)** Convert magnetometer measurements, roll, pitch, and magnetic declination into yaw. Input is magnetometer reading in microtesla and radians. Output is in radians in relation to True North.

```MATLAB
magut = [0.1, 0.1, 0.1];
roll = 0.1;
pitch = 0.1;
magvar_rad = -3.6;
psi = tilt_comp (magut, roll, pitch, magvar_rad);
```
