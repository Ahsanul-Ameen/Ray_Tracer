#include <bits/stdc++.h>
using namespace std;

#define SCENE_PATH "./scene.txt"

#define EPS	1e-9

//ambient, diffuse, specular, recursive reflection coefficient
#define AMB 	0
#define DIFF 	1
#define SPEC 	2
#define RRC 	3


#pragma once

extern int rec_level, dimensions;

class Object;
class Light;

extern vector <Object *> objects;
extern vector <Light *> lights;


ifstream scene_file(SCENE_PATH);




struct Point3D {
    double x, y, z;

    Point3D operator + (const Point3D& p) const {
        Point3D pp{};
        pp.x = x + p.x;
        pp.y = y + p.y;
        pp.z = z + p.z;
        return pp;
    }

    Point3D operator - (const Point3D& p) const {
        Point3D pp{};
        pp.x = x - p.x;
        pp.y = y - p.y;
        pp.z = z - p.z;
        return pp;
    }

    bool operator == (const Point3D& p) const {
        return x == p.x && y == p.y && z == p.z;
    }

    bool operator != (const Point3D& p) const {
        return !((*this) == p);
    }

    friend ostream &operator<<(ostream &out, const Point3D& p) {
		out << "(" << p.x << ", " << p.y << ", " << p.z << ")";
		return out;
	}

	friend istream &operator>>(istream &in, Point3D& p) {
		in >> p.x >> p.y >> p.z;
		return in;
	}
} eye;


struct Vector3D {
    Point3D endPoint;

    Vector3D() {}

	Vector3D(Point3D p) 
		: endPoint(p)
	{}  

    double getValue() const {
        return sqrt(endPoint.x * endPoint.x 
        			+ endPoint.y * endPoint.y 
        			+ endPoint.z * endPoint.z);
    }

    double dot(const Vector3D &v) const {
    	double dot_product = endPoint.x * v.endPoint.x
        				   		+ endPoint.y * v.endPoint.y
        				   		+ endPoint.z * v.endPoint.z;
    	return dot_product;
    }

    void normalize() {
        double val = getValue();
        endPoint.x /= val, endPoint.y /= val, endPoint.z /= val;
    }

    Vector3D operator * (double s) const {
        Vector3D v{};
        v.endPoint.x = endPoint.x * s;
        v.endPoint.y = endPoint.y * s;
        v.endPoint.z = endPoint.z * s;
        return v;
    }

    Vector3D operator * (const Vector3D& w) const {
        Vector3D v{};
        v.endPoint.x = endPoint.y * w.endPoint.z - endPoint.z * w.endPoint.y;
        v.endPoint.y = endPoint.z * w.endPoint.x - endPoint.x * w.endPoint.z;
        v.endPoint.z = endPoint.x * w.endPoint.y - endPoint.y * w.endPoint.x;
        return v;
    }

    Vector3D operator + (const Vector3D& v) const {
        Vector3D vv{};
        vv.endPoint = endPoint + v.endPoint;
        return vv;
    }

    Vector3D operator - (const Vector3D& v) const {
        Vector3D vv{};
        vv.endPoint = endPoint - v.endPoint;
        return vv;
    }

    friend ostream &operator<<(ostream &out, const Vector3D& p) {
		out << p.endPoint;
		return out;
	}
} u, r, l;


struct Color {
	double r, g, b; // stores value in range [0, 1]

	Color() 
		: r(0), g(0), b(0)
	{}

	Color(double r, double g, double b) 
		: r(r / 255.0), g(g / 255.0), b(b / 255.0)
	{
		clip();
	}

	Color(double *color)
		: r(color[0]), g(color[1]), b(color[2])
	{
		clip();
	}

	Color operator + (const Color &c) const {
        Color sum;
        sum.r = r + c.r;
        sum.g = g + c.g;
        sum.b = b + c.b;
    	sum.clip();
    	return sum;
    }

    Color operator * (const Color &c) const {
        Color prod;
        prod.r = r * c.r;
        prod.g = g * c.g;
        prod.b = b * c.b;
    	prod.clip();
    	return prod;
    }

    Color operator * (const double &coeff) const {
        Color scale;
        scale.r = r * coeff;
        scale.g = g * coeff;
        scale.b = b * coeff;
    	scale.clip();
    	return scale;
    }

    void clip() {
    	r = max(r, 0.0); r = min(r, 1.0);
    	g = max(g, 0.0); g = min(g, 1.0);
    	b = max(b, 0.0); b = min(b, 1.0);
    }

    void getColor(double *color) {
    	clip();
    	color[0] = r, color[1] = g, color[2] = b;
    }

};

class Ray {
public:
	Point3D start;
	Vector3D dir;
public:
	Ray(Point3D startPoint, Vector3D direction) 
		: start(startPoint), dir(direction)
	{
		dir.normalize();
	}

	Point3D getPoint(const double &t) const {
		return start + (dir * t).endPoint;
	}

	friend ostream &operator<<(ostream &out, const Ray& r) {
		out << "[ " << r.start << " ---> " << r.dir << " ]";
		return out;
	}
};



extern void drawSphere(double radius, int slices, int stacks);
extern void drawCheckerBoard(const Point3D &bottom_left, double tileWidth);
extern void drawTriangle(const vector<Point3D> &points);




class Light {
public:
	Point3D light_pos;
	double color[3];
public:
	Light(const Point3D &light_pos)
		: light_pos(light_pos)
	{}

	void setColor() {
		for(int i = 0; i < 3; i++)
			scene_file >> color[i];
	}

	void draw() {
		//cout << "Light :: draw()" << endl;
        glColor3f(color[0], color[1], color[2]);
        glPushMatrix(); {
            glTranslatef(light_pos.x, light_pos.y, light_pos.z);
            drawSphere(2, 15, 15);
        } glPopMatrix();
	}

	friend ostream &operator<<(ostream &out, const Light &l) {
		out << l.light_pos << endl;
		for( int i = 0; i < 3; i++ ) {
			out << l.color[i] << " ";
		}
		return out;
	}
};




class Object {
public:
	string type;
	vector<Point3D> reference_point; // should have x, y, z
	
	double height, width, length;
	double color[3];
	double coEfficients[4]; //ambient, diffuse, specular, recursive reflection coefficient
	// reflection coefficients
	int shine; // exponent term of specular component

public:
	virtual double getIntercectParam(Ray *r) = 0;

public:
	Object() {
		type = "object";
	}

	virtual void draw() = 0;	

	friend ostream &operator<<(ostream &out, const Object &o) {
		out << "Type : " << o.type << endl;
		out << "Color : ";
		for(int i = 0; i < 3; i++) out << o.color[i] << " ";
		out << "\nshine : " << o.shine << endl;
		return out;
	}

	virtual void setColor() {
		for(int i = 0; i < 3; i++)
			scene_file >> color[i];
	}

	virtual void setCoEfficients() {
		for(int i = 0; i < 4; i++)
			scene_file >> coEfficients[i];
	}

	virtual void setShine() {
		scene_file >> shine;
	}

	virtual Color getColorAt(const Point3D &intersectionPoint) {
		return Color();
	}

	virtual Vector3D getNormatAt(const Point3D &intersectionPoint) = 0;



	double intersect(Ray *r, double *color, int level) {

		// code for finding intersecting t min

		double t = getIntercectParam(r);



		if(level == 0)
			return t; 

		if (t < EPS) {
			cout << "ERROR t_min is almost -ve !" << endl;
		}

		Point3D intersectionPoint = r->getPoint(t);
		//cout << "{" << intersectionPoint << "}" << endl;

		Color intersectionPointColor = getColorAt(intersectionPoint);

		// color = intersectionPointColor * coEfficient[AMB]
		Color pointColor = intersectionPointColor * coEfficients[AMB];
		
		//calculate normal at intersectionPoint
		Vector3D normal = getNormatAt(intersectionPoint);
		normal.normalize();

		// MAKE SURE THAT THE NORMAL IS FACING THE RAY
		if(normal.dot(r->dir) > EPS) {
			normal = normal * (-1.0);
		}




		// for each light l in lights
		for(const auto &light : lights) {	
			// cast ray_l from l.light_pos to intersectionPoint
			
			Vector3D lray_dir({light->light_pos - intersectionPoint});
			double distance_from_source = lray_dir.getValue();

			// BE VERY CAREFUL ABOUT THE DIRECTION OF THIS TYPE OF RAY
			Ray lightRay(light->light_pos, lray_dir); // automatically normalized
			// also prevent internal reflections
			lightRay.start = lightRay.getPoint(0.47);


			bool isShadowed = false;

			for(const auto &o : objects) {
				double it = o->getIntercectParam(&lightRay);
				if (it > EPS && it + EPS < distance_from_source) { // bug fixed
					isShadowed = true;
					break;
				}
			}

			if(!isShadowed) {

                double lambertValue = max({lightRay.dir.dot(normal), 0.0});
                
                // R = L - 2(L.n)n 
                Vector3D reflected_ray = lightRay.dir - normal * (2.0 * lightRay.dir.dot(normal));	
                // oppose the direction to comput phong value
                reflected_ray = reflected_ray * (-1.0);

                reflected_ray.normalize();

                double phongValue = max({reflected_ray.dot(r->dir), 0.0});

                // l.color works as the source intensity, I_s here
                Color lightColor(light->color);

                // diffuse color 	color += l.color*coEfficient[DIFF]*lambertValue*intersectionPointColor
                pointColor = pointColor + ((lightColor * intersectionPointColor) * coEfficients[DIFF] * lambertValue); 

                // specular component	color += l.color*coEfficient[SPEC]*phongValue shine *intersectionPointColor
                pointColor = pointColor + ((lightColor * intersectionPointColor) * coEfficients[SPEC] * pow(phongValue, shine));
			
			}	
		}


		if (level < rec_level) {
			// construct reflected ray from intersection point
			// actually slightly forward from the point (by moving the
			// start a little bit towards the reflection direction)
			// to avoid self intersection

			// R = L - 2(L.n)n 		(BE CAREFUL HERE)
            Vector3D ref_camera_ray = r->dir - normal * (2.0 * r->dir.dot(normal)) ;
            ref_camera_ray.normalize();

            Ray reflectedCameraRay(intersectionPoint, ref_camera_ray);
            // prevent internal reflection
           	reflectedCameraRay.start = reflectedCameraRay.getPoint(0.47); 

            // find t_min from the nearest intersecting object, using
			// intersect() method, as done in the capture() method
			// if found, call intersect(r reflected , color reflected , level+1)
            {
				double *ref_color = new double[3];

				double nearest = -1;
				double t_min = INT_MAX, t;

				for(int k = 0, n = objects.size(); k < n; k++) {
					t = objects[k]->intersect(&reflectedCameraRay, ref_color, 0);
					if(t > EPS && t + EPS < t_min) {
						t_min = t;
						nearest = k;
					}		
				}

				if(nearest != -1) {
					t_min = objects[nearest]->intersect(&reflectedCameraRay, ref_color, level + 1);

					// color += color reflected * coEfficient[REC_REFFLECTION]	
		
					pointColor = pointColor + (Color(ref_color) * coEfficients[RRC]);				
				}

				delete[] ref_color;
				ref_color = nullptr;
			}
		}

		// -----------------------------

		pointColor.getColor(color); // produce the color vector and return

		return t;

	}
};


class Sphere : public Object {
public:
	double getIntercectParam(Ray *r) override {

		const Point3D &center = reference_point.back();
		const Point3D &startPoint = r->start;
		const Vector3D &dir = r->dir;
		const Vector3D &rc(startPoint - center);

		double a = dir.dot(dir);
		double b = 2.0 * dir.dot(rc);
		double c = rc.dot(rc) - (length * length);
		
		double d = b * b - 4.0 * a * c;
		
		if (d < 0) 
			return -1.0; 
		else {
			d = sqrt(d);
			double tp = (-b + d) / (2.0 * a);
			double tn = (-b - d) / (2.0 * a);

			if(tp >= 0 &&  tn >= 0) 
				return min(tp, tn);
			else 
				return max(tp, tn);
		}
	}

public:
	Sphere(Point3D center, double radius) {
		type = "sphere";
		reference_point.push_back(center);
		length = radius;
	}

	void draw() {
		//cout << "Sphere :: draw()" << endl;
		const Point3D p = reference_point.back();
        glColor3f(color[0], color[1], color[2]);
        glPushMatrix(); {
            glTranslatef(p.x, p.y, p.z);
            drawSphere(length, 30, 30);
        } glPopMatrix();
	}

	Color getColorAt(const Point3D &intersectionPoint) override {
		return Color(color);
	}

	Vector3D getNormatAt(const Point3D &intersectionPoint) override {
		Vector3D normal_vector(intersectionPoint - reference_point.back());
		normal_vector.normalize();
		return normal_vector;
	}
};


class Triangle : public Object {
public:
	double getIntercectParam(Ray *r) override {
		//cout << "triangle" << endl;

		const Point3D &v0 = reference_point[0];
		const Point3D &v1 = reference_point[1];
		const Point3D &v2 = reference_point[2];

		Vector3D edge_one(Point3D (v1 - v0));
		Vector3D edge_two(v2 - v0);
		Vector3D h = r->dir * edge_two;
		double a = edge_one.dot(h);
		if (a > -EPS && a < EPS) {
			return -1.0;	// This ray is parallel to this triangle.
		}

		double f = 1.0 / a;
		Vector3D s(r->start - v0);
		double u = f * s.dot(h);
		if (u < 0.0 || u > 1.0)
        	return -1.0;

        Vector3D q = s * edge_one;
        double v = f * r->dir.dot(q);
        if (v < 0.0 || u + v > 1.0)
        	return -1.0;
        
        // At this stage we can compute t to find out where the intersection point is on the line.
        double t = f * edge_two.dot(q);

        // (t <= EPS) means that there is a line intersection but not a ray intersection.
        return (t > EPS) ? t : -1.0;
	}

public:
	Triangle(Point3D a, Point3D b, Point3D c) {
		type = "triangle";
		reference_point.push_back(a);
		reference_point.push_back(b);
		reference_point.push_back(c);
	}

	void draw() {
		//cout << "Triangle :: draw()" << endl;
		glColor3f(color[0], color[1], color[2]);
        drawTriangle(reference_point); 
	}

	Color getColorAt(const Point3D &intersectionPoint) override {
		return Color(color);
	}

	Vector3D getNormatAt(const Point3D &intersectionPoint) override {
		const Point3D &v0 = reference_point[0];
		const Point3D &v1 = reference_point[1];
		const Point3D &v2 = reference_point[2];

		Vector3D edge_one(v1 - v0);
		Vector3D edge_two(v2 - v0);
		Vector3D normal_vector = edge_one * edge_two;
		normal_vector.normalize();

		return normal_vector;
	}

};

class General : public Object {
public:
	vector <double> quadratic_values;

public:

	bool isInsideBoundingBox(Ray *r, double t) {
		Point3D p1 = reference_point.back();
		Point3D p2 = p2 + Point3D({length, width, height});

		Point3D p = r->getPoint(t);

		if (length >= EPS) {
			if ( p.x < p1.x - EPS ) return false;
			if ( p.x > p2.x + EPS ) return false;
		}

		if (width >= EPS) {
			if ( p.y < p1.y - EPS ) return false;
			if ( p.y > p2.y + EPS ) return false;
		}

		if (height >= EPS) {
			if ( p.z < p1.z - EPS ) return false;
			if ( p.z > p2.z + EPS ) return false;
		}

	   	return true;
	}

	double getIntercectParam(Ray *r) override {
		//cout << "general" << endl;

		double A = quadratic_values[0];
		double B = quadratic_values[1];
		double C = quadratic_values[2];
		double D = quadratic_values[3];
		double E = quadratic_values[4];
		double F = quadratic_values[5];
		double G = quadratic_values[6];
		double H = quadratic_values[7];
		double I = quadratic_values[8];
		double J = quadratic_values[9];

		double rox = r->start.x, roy = r->start.y, roz = r->start.z;
		double rdx = r->dir.endPoint.x, rdy = r->dir.endPoint.y, rdz = r->dir.endPoint.z;


		const double a = A * rdx * rdx +
						 B * rdy * rdy +
						 C * rdz * rdz +
						 D * rdx * rdy + 
						 E * rdx * rdz + 
						 F * rdy * rdz;
		
		const double b = 2.0 * A * rox * rdx + 
						 2.0 * B * roy * rdy +
						 2.0 * C * roz * rdz + 
						 D * (rox * rdy + roy * rdx) +
						 E * (rox * rdz + roz * rdx) +
						 F * (roy * rdz + roz * rdy) + 
						 G * rdx + 
						 H * rdy + 
						 I * rdz;
			
		const double c = A * rox * rox + 
						 B * roy * roy +
						 C * roz * roz + 
						 D * rox * roy + 
						 E * rox * roz + 
						 F * roy * roz + 
						 G * rox + 
						 H * roy + 
						 I * roz + 
						 J;


		double d = b * b - 4 * a * c;

		if (abs(a) < EPS) { // a == 0
			double t = - c / b;
			bool in = isInsideBoundingBox(r, t);

			return (in) ? t : -1.0;
		} else {

			if (d < EPS) return -1.0; // no intersection
			else {
				d = sqrt(d);
				double tp = (-b + d) / (2.0 * a);
				double tn = (-b - d) / (2.0 * a);

				if(tp < tn) 
					swap(tp, tn);

				bool p_in = isInsideBoundingBox(r, tp);
				bool n_in = isInsideBoundingBox(r, tn);

				if (p_in && n_in) 
					return tn > 0 ? tn : tp;
				else if (p_in)
					return tp;
				else if (n_in)
					return tn;
				else return -1.0;
			}
		}
	}

public:
	General(vector<double> quadratic_values, Point3D crp, double l, double w, double h) {
		type = "general";
		this->quadratic_values = quadratic_values;
		reference_point.push_back(crp);
		length = l;
		width = w;
		height = h;
	}

	void draw() {}

	Color getColorAt(const Point3D &intersectionPoint) override {
		return Color(color);
	}

	Vector3D getNormatAt(const Point3D &intersectionPoint) override {
		double A = quadratic_values[0];
		double B = quadratic_values[1];
		double C = quadratic_values[2];
		double D = quadratic_values[3];
		double E = quadratic_values[4];
		double F = quadratic_values[5];
		double G = quadratic_values[6];
		double H = quadratic_values[7];
		double I = quadratic_values[8];
		double J = quadratic_values[9];

		double del_x = 2 * A * intersectionPoint.x
					 + D * intersectionPoint.y 
					 + E * intersectionPoint.z 
					 + G;

		double del_y = 2 * B * intersectionPoint.y 
					 + D * intersectionPoint.x
					 + F * intersectionPoint.z 
					 + H;

		double del_z = 2 * C * intersectionPoint.z 
					 + E * intersectionPoint.x
					 + F * intersectionPoint.y 
					 + I;
		Vector3D normal_vector({del_x, del_y, del_z});
		// Rn must be normalized and also we have to find the normal for surface facing the ray.
		normal_vector.normalize();

		return normal_vector;
	}

};


// declaration
class Floor : public Object {
public:
	double getIntercectParam(Ray *r) override {
		//cout << "floor" << endl;

		double numerator = (-r->start.z);
		double denominator = (r->dir.endPoint.z);

		double t = -1.0;
		if(abs(denominator) >= EPS) {
			t = numerator / denominator;
		} 

		return t >= EPS ? t : -1.0;
	}

public:
	Floor(double floorWidth, double tileWidth) {
		type = "floor";
		reference_point.push_back({-floorWidth / 2.0, -floorWidth / 2.0, 0});
		width = floorWidth;
		length = tileWidth; // store as length of tile
	}

	void setColor() override {
		color[0] = 1, color[1] = 1, color[2] = 1; // dummy body
	}

	void setCoEfficients() {
		coEfficients[AMB] = 0.5, coEfficients[DIFF] =  0.4, coEfficients[SPEC] = 0.7, coEfficients[RRC] = 0.3;
	}

	void setShine() {
		shine = 5;
	}

	void draw() {
		// write codes for drawing a checkerboard-like
		// floor with alternate colors on adjacent tiles
		drawCheckerBoard(reference_point.back(), length);
	}

	Color getColorAt(const Point3D &intersectionPoint) override {
		if(abs(intersectionPoint.x) > width / 2 || abs(intersectionPoint.y) > width / 2)
			return Color(0, 0, 0);
		int i = int(abs(intersectionPoint.x) / length) + (intersectionPoint.x < 0);
		int j = int(abs(intersectionPoint.y) / length) + (intersectionPoint.y < 0);
		int s = i + j;
		return (!(s & 1)) ? Color(0, 0, 0) : Color(255, 255, 255);
	}

	Vector3D getNormatAt(const Point3D &intersectionPoint) override {
		return Vector3D({0.0, 0.0, 1.0});
	}

};





extern vector <Object *> objects;
extern vector <Light *> lights;


void loadData(bool withFloor) {
	

	cout << "-----Reading from scene.txt-----" << endl;

	scene_file >> rec_level >> dimensions;

	int n_object;
	scene_file >> n_object;
	
	Object *temp = nullptr;
	for (int i = 0; i < n_object; ++i) {
		string object_name; 
		scene_file >> object_name;
		//cout << "object : " << i + 1 << " : " << object_name << endl;

		if (object_name == "sphere") {

			// read sphere data
			Point3D center;
			scene_file >> center;
			double radius;
			scene_file >> radius;
			
			temp = new Sphere(center, radius);
			temp->setColor();
			temp->setCoEfficients();
			temp->setShine();

			objects.push_back(temp);

		} else if (object_name == "triangle") {

			// read triangle data
			Point3D a, b, c;
			scene_file >> a >> b >> c;

			temp = new Triangle(a, b, c);
			temp->setColor();
			temp->setCoEfficients();
			temp->setShine();

			objects.push_back(temp);

		} else if (object_name == "general") {

			// read general data
			vector <double> quadratic_values(10);
			for(auto &i : quadratic_values)
				scene_file >> i;
			
			Point3D cube_referece_point;
			scene_file >> cube_referece_point;

			double length, width, height;
			scene_file >> length >> width >> height;

			temp = new General(quadratic_values, cube_referece_point, length, width, height);
			temp->setColor();
			temp->setCoEfficients();
			temp->setShine();

			objects.push_back(temp);

		} else cerr << "Error during reading " << object_name << " data" << endl;
		
	}
	
	if (withFloor) {

		// create a floor object
		temp = new Floor(1000, 20); // you can change this value
		temp->setColor();
		temp->setCoEfficients();
		temp->setShine();

		objects.push_back(temp);

	}
	temp = nullptr;


	int n_light;
	scene_file >> n_light;

	Light *light;
	for (int i = 0; i < n_light; ++i) {
		Point3D l_pos;
		scene_file >> l_pos;

		light = new Light(l_pos);
		light->setColor();

		lights.push_back(light);
	}
	light = nullptr;

	cout << "-----Reading complete-----" << endl;
}



