
import numpy as np
from types import SimpleNamespace

NAD83_WGS84_dict = {
    "a":6378137,
    "b":6356752.3142,
    "flat":1/298.257223563,
    "a_dscp":"Equatorial Radius, meters",
    "b_dscp":"Polar Radius, meters",
    "flat_dscp":"Flattening (a-b)/a",
}

NAD83_WGS84_obj = SimpleNamespace(**NAD83_WGS84_dict)

class BasicMapObj:

    def __init__(self,ellipsoid=NAD83_WGS84_obj,zone=None):

        b=ellipsoid.b
        a=ellipsoid.a
        e=np.sqrt(1-b**2/a**2)

        self.b=b
        self.a=a
        self.e=e
        self.ep2=(e*a/b)**2
        self.n=(a-b)/(a+b)
        self.k0=0.9996
        self.easting = 500000
        self.zone = zone
        #self.northing = None

        if zone is None:
            self.lon0=None
            self.northing = None
        else:
            if zone.startswith('326'):
                zone = zone[3:5] + 'N'
                self.zone = zone
            elif zone.startswith('327'):
                zone = zone[3:5] + 'S'
                self.zone = zone

            if str(zone)[-1:].isnumeric(): # default is N, not S
                zone_number = int(zone)
                self.northing = 0
            else:
                zone_number = int(zone[:-1])
                if zone[-1] in ('N','n'):
                    self.northing = 0
                elif zone[-1] in ('S','s'):
                    self.northing = 1e7

            self.lon0 = (zone_number - 1)*6 -180 +3 # in Degrees  

    def calc_rho(self,lat_rad):
        a=self.a
        #b=self.b
        e=self.e

        return a*(1-e**2)/((1-e**2*(np.sin(lat_rad))**2)**(3/2))

    def calc_nu(self,lat_rad):
        a=self.a
        e=self.e
        return a / (1-(e*np.sin(lat_rad))**2)**0.5

    def calc_p(self,lon_rad):
        return lon_rad - np.radians(self.lon0)

    def calc_S(self,lat_rad):
        #S is the meridional arc
        a=self.a
        n=self.n
        a_p = 1  *    a *          (1 - n + 5/4*(n**2-n**3) + 81/64*(n**4-n**5))
        b_p = 3/2 *   a *    n   * (1 - n + 7/8*(n**2-n**3) + 55/64*(n**4))
        c_p = 15/16 * a * (n**2) * (1 - n + 3/4*(n**2-n**3))
        d_p = 35/48 * a * (n**3) * (1 - n + 11/16*(n**2))
        e_p = 315/512*a * (n**4) * (1 - n)

        s =  a_p*lat_rad  \
           - b_p*np.sin(2*lat_rad)  \
           + c_p*np.sin(4*lat_rad) \
           - d_p*np.sin(6*lat_rad) \
           + e_p*np.sin(8*lat_rad) \

        return s

    def calc_K3(self,nu,lat_rad):
        k0 = self.k0
        ep2 = self.ep2

        k_3 = k0*nu*np.sin(lat_rad)* (np.cos(lat_rad))**3 / 24
        k_3 *= 5 - (np.tan(lat_rad))**2 + 9 * ep2 * (np.cos(lat_rad))**2 + 4 * (ep2**2) * (np.cos(lat_rad))**4

        return k_3

    def calc_K5(self,nu,lat_rad):
        k0 = self.k0
        ep2 = self.ep2

        k_5 = k0 * nu * (np.cos(lat_rad))**3 /6
        k_5 *= 1 - (np.tan(lat_rad))**2 + ep2 * (np.cos(lat_rad))**2

        return k_5

    def estimate_lon0(self, lon_deg):
        if self.lon0 is None:
            major_lon = np.median(lon_deg)
            central_meridians = np.arange(0,60,1)*6 - 180 +3
            close_meridian = central_meridians[np.argmin(np.abs(major_lon-central_meridians))]
            self.lon0 = close_meridian
            self.zone =  int((close_meridian-3 +180)/6)+1  #(zone_number - 1)*6 -180 +3
        else:
            #use  lon0 during initialization
            pass

    def estimate_northing(self,lat_deg):
        if self.northing is None:
            major_lat = np.median(lat_deg)
            if major_lat>0:
                self.northing=0
            else:
                self.northing=1e7

    def convert_xycoord(self,lat_deg,lon_deg):
        lat_rad = np.radians(lat_deg)
        lon_rad = np.radians(lon_deg)

        self.estimate_lon0(lon_deg)
        #print(self.lon0)

        self.estimate_northing(lat_deg)

        s = self.calc_S(lat_rad)
        k0 = self.k0
        nu = self.calc_nu(lat_rad)
        p = self.calc_p(lon_rad)

        k_1 = s*k0
        k_2 = k0*nu*np.sin(2*lat_rad)/4
        k_3 = self.calc_K3(nu,lat_rad)

        y = k_1 + k_2 * (p**2) + k_3 * (p**4) + self.northing

        k_5 = self.calc_K5(nu,lat_rad)
        k_4 = k0 * nu * np.cos(lat_rad)

        x = k_4*p + k_5*(p**3)+ self.easting

        return x,y

    ########################
    #https://gdal.org/en/stable/proj_list/transverse_mercator.html
    # ref: Snyder J.P. (1987) Map projections a working manual, U.S. Geological Survey Professional Paper 1395, 1987. page.61
    def convert_xycoord_gdal(self, lat_deg,lon_deg):
        lat_rad = np.radians(lat_deg)
        lon_rad = np.radians(lon_deg)

        self.estimate_lon0(lon_deg)
        self.estimate_northing(lat_deg)

        k0 = self.k0
        E = (self.e)**2
        p = self.calc_p(lon_rad)
        cos_lat = np.cos(lat_rad)
        sin_lat = np.sin(lat_rad)
        tan_lat = sin_lat / cos_lat
        tan2_lat = tan_lat**2

        e_p2 = self.ep2

        nu = self.calc_nu(lat_rad)
        #nu = self.a / np.sqrt(1 - E * sin_lat**2)
        C = e_p2 * cos_lat**2
        A = cos_lat * p

        E2=E**2
        E3=E**3

        M1 = 1 - E / 4 - 3 * E2 / 64 - 5 * E3 / 256
        M2 = 3 * E / 8 + 3 * E2 / 32 + 45 * E3 / 1024
        M3 = 15 * E2 / 256 + 45 * E3 / 1024
        M4 = 35 * E3 / 3072
        M = self.a * (M1 * lat_rad -
             M2 * np.sin(2 * lat_rad) +
             M3 * np.sin(4 * lat_rad) -
             M4 * np.sin(6 * lat_rad))

        #M = a[(1 - e2/4 - 3e4/64 - 5e6/256 -....)* - (3e2/8 + 3e4/32 + 45e6/1024+....)sin2*
		#+ (15e4/256 + 45e6/1024 +.....)sin4* - (35e6/3072 + ....)sin6* + .....]

        x = k0 * nu * (A +
                            A**3 / 6 * (1 - tan2_lat + C) +
                            A**5 / 120 * (5 - 18 * tan2_lat + tan2_lat**2 + 72 * C - 58 * e_p2))+ self.easting

        y = k0 * (M + nu * tan_lat * (A**2 / 2 +
                                        A**4 / 24 * (5 - tan2_lat + 9 * C + 4 * C**2) +
                                        A**6 / 720 * (61 - 58 * tan2_lat + tan2_lat**2 + 600 * C - 330 * e_p2)))+ self.northing

        return x,y

    ########################

    def calc_mu(self):  #calc_e1_mu(self):
        e=self.e
        a=self.a

        mu_recip = a * (1-0.25*(e**2) -3/64*(e**4) -5/256 * (e**6)) 
        #e1 = (1 - eee) / (1 + eee) # same as self.n
        return mu_recip

    # ref : Snyder J.P. (1987) Map projections a working manual, U.S. Geological Survey Professional Paper 1395, 1987.  page.63
    # https://pubs.usgs.gov/pp/1395/report.pdf
    def convert_latlon(self,x,y):
        x_in = x - self.easting
        y_in = y - self.northing

        ep2 = self.ep2
        a = self.a
        e =self.e

        k0 = self.k0

        M = y_in / k0

        mu_recip =  self.calc_mu()  #self.calc_e1_mu()
        e1=self.n
        mu = M / mu_recip

        J1 = 3/2 * e1 - 27/32 * (e1**3)
        J2 = 21/16*(e1**2) -55/32*(e1**4)
        J3 = 151/96 *   (e1**3)
        J4 = 1097/512 * (e1**4)

        fp = mu + J1*np.sin(2*mu) + J2*np.sin(4*mu) + J3*np.sin(6*mu) + J4*np.sin(8*mu)

        C1 = ep2*(np.cos(fp))**2
        T1 = (np.tan(fp))**2
        R1 = a*(1-e**2) / (1-(e*np.sin(fp))**2)**1.5
        N1 = a / (1-(e*np.sin(fp))**2)**0.5
        D = x_in / N1 / k0

        Q1 = N1*np.tan(fp)/R1
        Q2 = D**2 / 2
        Q3 = (5 + 3*T1 + 10*C1 - 4*C1**2 -9*ep2) * D**4 / 24
        Q4 = (61 + 90*T1 + 298*C1 +45*T1**2 - 3*C1**2 -252*ep2) * D**6 /720

        lat_out = fp - Q1*(Q2-Q3+Q4)

        Q5 = D
        Q6 = (1 + 2*T1 + C1) * D**3 / 6
        Q7 = (5 - 2*C1 + 28*T1 -3*C1**2 + 8*ep2 +24*T1**2) * D**5 / 120

        lon_out = np.radians(self.lon0) + (Q5-Q6+Q7) / np.cos(fp)

        return np.degrees(lat_out), np.degrees(lon_out)
