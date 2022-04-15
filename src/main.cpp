/* \author Aaron Brown */
// Create simple 3d highway enviroment using PCL
// for exploring self-driving car sensors

#include "highway.h"
#include "gnuplot-iostream.h"

int main(int argc, char** argv)
{
    Gnuplot gp;
	pcl::visualization::PCLVisualizer::Ptr viewer(new pcl::visualization::PCLVisualizer("3D Viewer"));
	viewer->setBackgroundColor(0, 0, 0);

	// set camera position and angle
	viewer->initCameraParameters();
	float x_pos = 0;
	viewer->setCameraPosition ( x_pos-26, 0, 15.0, x_pos+25, 0, 0, 0, 0, 1);

	Highway highway(viewer);

	//initHighway(viewer);

	int frame_per_sec = 30;
	int sec_interval = 10;
	int frame_count = 0;
	int time_us = 0;

	double egoVelocity = 25;

	while (frame_count < (frame_per_sec*sec_interval))
	{
		viewer->removeAllPointClouds();
		viewer->removeAllShapes();

		//stepHighway(egoVelocity,time_us, frame_per_sec, viewer);
		highway.stepHighway(egoVelocity,time_us, frame_per_sec, viewer);
		viewer->spinOnce(1000/frame_per_sec);
		frame_count++;
		time_us = 1000000*frame_count/frame_per_sec;
		
	}

    //TODO must optimize plotting

    float radar_thold = 7.815;
    float lidar_thold = 5.991;

    // Iterate over car in traffic
    for(auto & elem : highway.traffic)
    {
        if (elem.ukf.nis_check_) {
            std::vector<std::pair<double, double> > xy_pts_radar_nis;
            // std::vector<std::pair<double, double> > xy_pts_lidar_nis;
            for(double i = 0; i < elem.ukf.nis_radar_.size(); ++i) {
                xy_pts_radar_nis.push_back(std::make_pair(i, elem.ukf.nis_radar_[i]));
                // xy_pts_lidar_nis.push_back(std::make_pair(i, elem.ukf.nis_lidar_[i]));
            }

            std::vector<std::pair<double, double> > xy_pts_radar_thold;
            xy_pts_radar_thold.push_back(std::make_pair(0, radar_thold));
            xy_pts_radar_thold.push_back(std::make_pair(elem.ukf.nis_radar_.size(), radar_thold));

            gp << "plot" << gp.file1d(xy_pts_radar_nis) << "with lines title 'Radar NIS',"
               << gp.file1d(xy_pts_radar_thold) << "with lines title 'Radar NIS threshold'" << std::endl;

            // std::vector<std::pair<double, double> > xy_pts_lidar_thold;
            // xy_pts_lidar_thold.push_back(std::make_pair(0, lidar_thold));
            // xy_pts_lidar_thold.push_back(std::make_pair(elem.ukf.nis_lidar_.size(), lidar_thold));

            // gp << "plot" << gp.file1d(xy_pts_lidar_nis) << "with lines title 'Lidar NIS',"
               // << gp.file1d(xy_pts_lidar_thold) << "with lines title 'Lidar NIS threshold'" << std::endl;
        }
    }

}