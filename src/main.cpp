#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"

// for spline fitting:
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "spline.h"

using namespace std;
using namespace Eigen;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Global starting lane
int lane = 1; // 0: left, 1: center, 2: right
int old_lane = 1; // to remember predecesor lane while changing lanes

// Global reference velocity
double ref_vel = 0.0; // mph

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}

int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
{
	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}
	}
	return closestWaypoint;
}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}


int transition_function(double car_s, int prev_size, vector<vector<double>> sensor_fusion){
  // Try to introduce lane change using a cost function
  // Only consider states which can be reached from current FSM state.
  vector<int> next_lanes;
  int target_lane;

  if(lane == 0){
    next_lanes.push_back(1);
  }
  if(lane == 1){
    next_lanes.push_back(0);
    next_lanes.push_back(2);
  }
  if(lane == 2){
    next_lanes.push_back(1);
  }

  // Keep track of the total cost of each state.
  vector<int> cost;
  cost.push_back(999);
  cost.push_back(999);
  cost.push_back(999);

  for(int l = 0; l < next_lanes.size(); l++){
    target_lane = next_lanes[l];
    for(int i = 0; i < sensor_fusion.size(); i++){
      // Car is in next_lanes / possible successor lane
      float d = sensor_fusion[i][6];
      double check_car_s = sensor_fusion[i][5];
      if(d < (2+4*target_lane+2) && d > (2+4*target_lane-2) && (fabs(check_car_s-car_s) > 60)){
        double vx = sensor_fusion[i][3];
        double vy = sensor_fusion[i][4];
        double check_speed = sqrt(vx*vx+vy*vy);

        // Using previous points and velocity to predict position in the future
        check_car_s += ((double)prev_size*0.02*check_speed);

        // Lane Change Cost Functions:
        // Car is ahead:
        if(check_car_s >= car_s){
          // Car is far away
          if((fabs(check_car_s-car_s)) > 50){
            cost[target_lane] -= 500;
          }
          // Car is too close
          else{
            cost[target_lane] -= 0;
          }
        }
        // Car is behind:
        else{
          // Car is far away
          if((fabs(check_car_s-car_s)) > 30){
            if(check_speed < 15){
              cost[target_lane] -= 500;
            }
            // Cars is not too close but too fast
            else{
              cost[target_lane] -= 0;
            }
          }
          // Car is too close
          else{
            cost[target_lane] -= 0;
          }
        }
      }
    }
  }

  // Find minimum cost state
  int best_next_state = lane;
  int min_cost = 999;
  for(int l = 0; l < next_lanes.size(); l++){
    target_lane = next_lanes[l];
    int target_cost = cost[target_lane];
    if(target_cost < min_cost){
      min_cost = target_cost;
      best_next_state = target_lane;
    }
  }

  // Choose best next state / lane
  return best_next_state;
}

int main() {
  uWS::Hub h;

  // Path Planner is initialized here!
  // PathPlaner pathplaner;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	// Note: Returns the previous list but with processed points removed, can be a nice tool to show how far along the path has processed since last time.
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	// auto sensor_fusion = j[1]["sensor_fusion"];
          	// use double instead of auto in order to be able to transfer to transition function
          	vector<vector<double>> sensor_fusion = j[1]["sensor_fusion"];

          	json msgJson;

          	vector<double> next_x_vals;
          	vector<double> next_y_vals;

          	int prev_size = previous_path_x.size();

          	// LANE CHANGE & COLLISSION AVOIDANCE:
          	// Use sensor fusion data of close by vehicles for deciding on lane change and velocity of ego car

          	bool change_lane = false;
          	bool slow_down = false;
          	double sync_speed;
          	float d;
          	double vx;
          	double vy;
          	double check_speed;
          	double check_car_s;
          	int selected_lane;

            if(prev_size > 0){
              car_s = end_path_s;
            }

          	// Detect slower car ahead and trigger reaction:

          	for(int i = 0; i < sensor_fusion.size(); i++){
          	  d = sensor_fusion[i][6];
          	  // Check if car is in ego vehicle's lane
          	  if(d < (2+4*lane+2) && d > (2+4*lane-2) ){
          	    vx = sensor_fusion[i][3];
          	    vy = sensor_fusion[i][4];
                check_car_s = sensor_fusion[i][5];
                check_speed = sqrt(vx*vx+vy*vy);

          	    // Using previous points and velocity to predict position in the future
          	    check_car_s += ((double)prev_size*0.02*check_speed);
          	    // Check if car is ahead (>s value) and gap is VERY small
                if((check_car_s > car_s) && ((check_car_s-car_s) < 25)){
                  // If car in front is too close focus on increasing distance before considering lane change
                  slow_down = true;
                  // Save speed of slower car in the front of ego vehicle for later velocity adjustments
                  sync_speed = check_speed;
                }
                // Check if car is ahead (>s value) and too close
          	    else if((check_car_s > car_s) && ((check_car_s-car_s) < 40)){
          	      // Trigger changing lanes
          	      change_lane = true;
          	    }
          	  }
          	}

          	// Introduce Lane Change:

          	if (change_lane){
          	  // Use transition function (FSM) to decide for behavior / lane to change to
              selected_lane = transition_function(car_s, prev_size, sensor_fusion);
          	  // Safety check: is selected lane safe to change to?
              for(int i = 0; i < sensor_fusion.size(); i++){
                d = sensor_fusion[i][6];
                // Check if car is in selected lane
                if(d < (2+4*selected_lane+2) && d > (2+4*selected_lane-2) ){
                  vx = sensor_fusion[i][3];
                  vy = sensor_fusion[i][4];
                  check_car_s = sensor_fusion[i][5];
                  check_speed = sqrt(vx*vx+vy*vy);

                  // Using previous points and velocity to predict position in the future
                  check_car_s += ((double)prev_size*0.02*check_speed);
                  // Car is closer than 20 m
                  if(fabs(check_car_s-car_s) < 20){
                    // Stay in old lane, new lane is not safe!
                    selected_lane = lane;
                  }
                  // Car is 20 to 100 m close and faster!
                  if((fabs(check_car_s-car_s) > 20) && (check_speed > car_speed) && (fabs(check_car_s-car_s) < 100)){
                    // Stay in old lane, avoid annoying faster cars on the Autobahn!
                    selected_lane = lane;
                  }
                }
              }
              // Lane safe and approved for lane change
              lane = selected_lane;
          	}

          	// Status changing lanes:

          	if (lane != old_lane){
              // Additional Safety check to prevent collision during lane change
              for(int i = 0; i < sensor_fusion.size(); i++){
                d = sensor_fusion[i][6];
                // Check if car is in the same lane in which ego vehicle is changing into
                if(d < (2+4*lane+2) && d > (2+4*lane-2) ){
                  vx = sensor_fusion[i][3];
                  vy = sensor_fusion[i][4];
                  check_car_s = sensor_fusion[i][5];
                  check_speed = sqrt(vx*vx+vy*vy);

                  // Using previous points and velocity to predict position in the future
                  check_car_s += ((double)prev_size*0.02*check_speed);
                  // Car is less than 50 m behind ego vehicle and faster
                  if(((check_car_s-car_s) < -50) && (check_speed > car_speed)){
                    // Go back to old lane to avoid collision
                    lane = old_lane;
                  }
                  // Check if lane change is completed
                  if(car_d < (2+4*lane+2) && car_d > (2+4*lane-2) ){
                    // Forget the old lane
                    old_lane = lane;
                  }
                }
              }
          	}

          	// Velocity control:

          	// Initialize speed adjustments if car ahead is very close
          	if(slow_down){
          	  // Check if ego vehicle is faster than slower car ahead
          	  if (ref_vel > sync_speed){
                // Slow down approx. 5 m/s^2
          	    ref_vel -= 0.224;
          	  }
          	  // Check if car is very close but ego vehicle is slower (i.e. after lane change)
          	  else if (ref_vel < sync_speed){
          	    // Speed up slowly
          	    ref_vel += 0.1;
          	  }
          	}
          	// Check if there is no car ahead very close and ego vehicle is slower than reference velocity
          	else if (ref_vel < 49.5){
          	  // Speed up approx. 5 m/s^2
          	  ref_vel += 0.224;
          	}

          	// KEEPING THE LANE:
          	// Create a list of waypoints (cartesian) with a distance of 30 m
          	// In order to control velocity a spline function is introduced afterwards to interpolate additional points between those anchor points

          	// Get car's position (or previous path's end point) as starting reference
          	double ref_x = car_x;
          	double ref_y = car_y;
          	double ref_yaw = deg2rad(car_yaw);

            vector<double> ptsx;
            vector<double> ptsy;

          	// If previous path is almost empty, use car's position as starting reference
          	if(prev_size < 2){
          	  // Use points that make the path tangent to the car
          	  double prev_car_x = car_x - cos(car_yaw);
              double prev_car_y = car_y - sin(car_yaw);

              ptsx.push_back(prev_car_x);
              ptsx.push_back(car_x);
              ptsy.push_back(prev_car_y);
              ptsy.push_back(car_y);
          	}
          	// If previous path is already recorded, use its end point as starting reference
          	else{
          	  // Set previous path's end point as reference state
          	  ref_x = previous_path_x[prev_size-1];
          	  ref_y = previous_path_y[prev_size-1];

          	  double ref_x_prev = previous_path_x[prev_size-2];
          	  double ref_y_prev = previous_path_y[prev_size-2];
          	  ref_yaw = atan2(ref_y-ref_y_prev, ref_x-ref_x_prev);

          	  // Use points that make the path tangent to previous path's end point
          	  ptsx.push_back(ref_x_prev);
          	  ptsx.push_back(ref_x);
              ptsy.push_back(ref_y_prev);
              ptsy.push_back(ref_y);
          	}

          	// Add evenly 30 m spaced points ahead of the starting reference (Frenet)
          	vector<double> next_anchor_point_0 = getXY(car_s+30, (2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
            vector<double> next_anchor_point_1 = getXY(car_s+60, (2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
            vector<double> next_anchor_point_2 = getXY(car_s+90, (2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);

            ptsx.push_back(next_anchor_point_0[0]);
            ptsx.push_back(next_anchor_point_1[0]);
            ptsx.push_back(next_anchor_point_2[0]);

            ptsy.push_back(next_anchor_point_0[1]);
            ptsy.push_back(next_anchor_point_1[1]);
            ptsy.push_back(next_anchor_point_2[1]);

            for(int i = 0; i < ptsx.size(); i++){
              // Shift car's reference angle to 0 degrees
              double shift_x = ptsx[i]-ref_x;
              double shift_y = ptsy[i]-ref_y;

              ptsx[i] = shift_x*cos(0-ref_yaw) - shift_y*sin(0-ref_yaw);
              ptsy[i] = shift_x*sin(0-ref_yaw) + shift_y*cos(0-ref_yaw);
            }

            // Create a spline
            tk::spline s;

            // Set list of points (cartesian) to the spline
            s.set_points(ptsx, ptsy);

            // Begin with previous path points
            for(int i = 0; i < prev_size; i++){
              next_x_vals.push_back(previous_path_x[i]);
              next_y_vals.push_back(previous_path_y[i]);
            }

            // In order to achieve target reference velocity spline point have to be broken up
            double target_x = 30.0;
            double target_y = s(target_x);
            // Euclidean distance from 0
            double target_dist = sqrt((target_x)*(target_x) + (target_y)*(target_y));
            double x_add_on = 0;

            // Fill up the rest of the trajectory after filling it with previous points
            // Goal: Output 50 points
            for(int i = 1; i <= 50-prev_size; i++){
              // Car moves 50 times a sec to the next waypoint -> 20 ms
              double N = target_dist/(0.02*ref_vel/2.236936);
              double x_point = x_add_on+(target_x)/N;
              double y_point = s(x_point);

              x_add_on = x_point;

              double x_ref = x_point;
              double y_ref = y_point;

              // Rotate back to normal after rotating it earlier
              x_point = x_ref*cos(ref_yaw) - y_ref*sin(ref_yaw);
              y_point = x_ref*sin(ref_yaw) + y_ref*cos(ref_yaw);

              x_point += ref_x;
              y_point += ref_y;

              next_x_vals.push_back(x_point);
              next_y_vals.push_back(y_point);
            }

            // Output to simulator:
          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
















































































