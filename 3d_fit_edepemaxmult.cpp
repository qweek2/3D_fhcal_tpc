#include <iostream>
#include <fstream>
#include <vector>
#include <TFile.h>
#include <TNtuple.h>
#include <TLine.h>
#include <TGraph.h>
#include <TF1.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include "TCutG.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <stdio.h>
// #include "utils.h"
#include "TString.h"
#include "TH1D.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TH2.h"
// Load the ROOT libraries
#include "TH2F.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TLine.h"
#include <TVector3.h>
#include <TH3F.h>

vector<double> read_lines(int n, const string &filename)
{
    // Open the file for reading
    ifstream file(filename);
    // Check if the file was opened successfully
    if (!file.is_open())
    {
        cout << "Error opening file!" << endl;
        return {};
    }
    // Declare a vector to store the values
    vector<double> values;
    // Read the lines from the file
    double value1, value2, value3;
    for (int i = 1; i <= n + 1; i++)
    {
        if (file >> value1 >> value2 >> value3)
        {
            if (i == n || i == n + 1)
            {
                values.push_back(value1);
                values.push_back(value2);
                values.push_back(value3);
            }
        }
        else
        {
            cout << "Error reading file!" << endl;
            return {};
        }
    }
    // Close the file
    file.close();
    // Return the values
    return values;
}

vector<double> findMatchingLines(int i, int j, int k, int h_counter, string filename)
{
    vector<double> result;
    ifstream file(filename);
    ofstream output("output.txt", ios::app); // Open output file in append mode

    if (file.is_open())
    {
        int a, b, c;
        double d;
        int lineNum = 0; // Line number counter
        while (file >> a >> b >> c >> d)
        {
            if (a == i & b == j & c == k)
            {
                result.push_back(d);                                                                               // b or npart
                output << a << " " << b << " " << c << " " << d << " " << lineNum << " " << h_counter - 1 << endl; // Write to output file
            }
            lineNum++;
        }
        file.close();
        output.close(); // Close output file
    }
    else
    {
        cout << "Error opening file: " << filename << endl;
    }
    return result;
}

void fit()
{
    ofstream myfile, myfile2;
    myfile.open("3d_fit_list.txt");
    myfile2.open("log_3D.txt");

    // Open the ROOT file containing the TH3F histogram
    TFile *file = new TFile("mult_edep_test.root");
    if (!file || file->IsZombie())
    {
        std::cerr << "Failed to open the ROOT file." << std::endl;
    }

    // Access the TH3F histogram from the ROOT file
    TH3F *histogram = dynamic_cast<TH3F *>(file->Get("EdepEmaxMult"));
    if (!histogram)
    {
        std::cerr << "Failed to retrieve the TH3F histogram from the ROOT file." << std::endl;
    }
    cout << "INTEGRAL:  " << histogram->Integral() << endl;
    TH3F *EdepEmaxMult_1 = new TH3F("EdepEmaxMult_1", "EdepEmaxMult_1", 250, 0, 30, 200, 0, 5, 300, 0, 300);

    TH3F *EdepEmaxMult[20];
    for (int i = 0; i < 20; i++)
    {
        EdepEmaxMult[i] = new TH3F(TString::Format("EdepEmaxMult%d", i), TString::Format("EdepEmaxMult%d", i), 250, 0, 30, 200, 0, 5, 300, 0, 300);
    }

    TH3F *tmp_for_fit3 = new TH3F("tmp_for_fit3", "tmp_for_fit3", 250, 0, 30, 200, 0, 5, 300, 0, 300);

    TH1F *ImpPar[24];
    for (int i_impact = 0; i_impact < 25; i_impact++)
    {
        ImpPar[i_impact] = new TH1F(TString::Format("impact_%d", i_impact), TString::Format("impact_%d", i_impact), 500, 0., 16.);
    }

    TH1F *h_npart[24];
    for (int i = 0; i < 25; i++)
    {
        h_npart[i] = new TH1F(TString::Format("npart_%d", i), TString::Format("npart_%d", i), 500, 0., 500.);
    }
    // ellipse parameters
    double a_ell = 9.1657;
    double b_ell = 0.7625;
    double th = 0.148095;
    double x0 = 5.1264;
    double y0 = 0.1125;

    double ys_d, ys_rot_d, xs_rot_d, ys_rot, ys, xs_rot, zm;
    for (double xs_d = -20; xs_d < 40; xs_d = xs_d + 0.1)
    {
        double xs = xs_d;
        ys = (b_ell * sqrt(a_ell * a_ell - xs * xs + 2 * xs * x0 - x0 * x0)) / a_ell + y0;
        xs_rot = x0 + (xs - x0) * cos(th) - (ys - y0) * sin(th); // получили новые коорд в повернутой СК
        ys_rot = y0 + (xs - x0) * sin(th) + (ys - y0) * cos(th);

        double delta = 0, delta_prev = 100000;
        double value1, value2;
        ifstream file("3d_fit_list_pol_up.txt");
        if (file.is_open())
        {
            if (xs_rot < 0.879665)
                continue;
            while (file >> value1 >> value2)
            {
                if (xs_rot > 0 & ys_rot > 0)
                {
                    delta = fabs(value2 - xs_rot);
                    // cout << "x " << xs_rot << " val txt " << value2 << " delta prev " << delta_prev << " delta now " << delta << endl;
                    if (delta > delta_prev)
                    {
                        zm = value1;
                        cout << xs_rot << " " << ys_rot << " " << zm << endl;
                        tmp_for_fit3->Fill(xs_rot, ys_rot, zm, 1); // up
                        break;
                    }
                    else
                        delta_prev = delta;
                }
            }
            file.close();
        }
        else
        {
            cout << "Unable to open file." << endl;
        }

        if (xs_rot > 0 & ys_rot > 0)
        {
            // polin4 = 0.879665 + ys_rot * 0.362493 + ys_rot * ys_rot * (-0.00369157) + ys_rot * ys_rot * ys_rot * 1.57008e-05 + ys_rot * ys_rot * ys_rot * ys_rot * -2.49992e-08;
            myfile << "up " << xs_rot << " " << ys_rot << endl;
        }
    }
    // down
    for (double xs_d = -20; xs_d < 40; xs_d = xs_d + 0.1)
    {
        ys_d = y0 - (b_ell * sqrt(a_ell * a_ell - xs_d * xs_d + 2 * xs_d * x0 - x0 * x0)) / a_ell;
        xs_rot_d = x0 + (xs_d - x0) * cos(th) - (ys_d - y0) * sin(th); // получили новые коорд в повернутой СК
        ys_rot_d = y0 + (xs_d - x0) * sin(th) + (ys_d - y0) * cos(th);

        double delta = 0, delta_prev = 100000;
        double value1, value2;
        ifstream file("3d_fit_list_pol_down.txt");
        if (file.is_open())
        {
            while (file >> value1 >> value2)
            {
                if (xs_rot_d > 0 & ys_rot_d > 0)
                {
                    delta = fabs(value2 - xs_rot_d);
                    // cout << "x " << xs_rot << " val txt " << value2 << " delta prev " << delta_prev << " delta now " << delta << endl;
                    if (delta > delta_prev)
                    {
                        zm = value1;
                        cout << xs_rot_d << " " << ys_rot_d << " " << zm << endl;
                        tmp_for_fit3->Fill(xs_rot_d, ys_rot_d, zm, 1); // down
                        break;
                    }
                    else
                        delta_prev = delta;
                }
            }
            file.close();
        }
        if (xs_rot_d > 0 & ys_rot_d > 0)
        {
            // polin4 = 0.879665 + ys_rot_d * 0.362493 + ys_rot_d * ys_rot_d * (-0.00369157) + ys_rot_d * ys_rot_d * ys_rot_d * 1.57008e-05 + ys_rot_d * ys_rot_d * ys_rot_d * ys_rot_d * -2.49992e-08;
            myfile << "down " << xs_rot_d << " " << ys_rot_d << endl;
        }
    }
    /*for (double xs_d = 0; xs_d < 300; xs_d = xs_d + 0.01)
    {
        polin4 = 0.879665 + xs_d * 0.362493 + xs_d * xs_d * (-0.00369157) + xs_d * xs_d * xs_d * 1.57008e-05 + xs_d * xs_d * xs_d * xs_d * -2.49992e-08;
        // myfile2 << xs_d << " " << polin4 << endl;
    }*/

    int binCount = 0;
    int h_counter = 0;
    for (int line = 1; line < 144; line++)
    {
        vector<double> values = read_lines(line, "points_3D.txt");
        if (values.size() == 6)
        {
            cout << endl;
            cout << "Line " << line << ": " << values[0] << " " << values[1] << " " << values[2] << endl;
            cout << "Line " << line + 1 << ": " << values[3] << " " << values[4] << " " << values[5] << endl;
        }

        // Define the coordinates of the two points
        double x1 = values[0];
        double y1 = values[1];
        double z1 = values[2];
        double x2 = values[3];
        double y2 = values[4];
        double z2 = values[5];

        // 3.50577 0.630896 7.86
        // 3.60447 0.646994 8.18

        //  Calculate the direction vector of the line
        TVector3 lineVector(x2 - x1, y2 - y1, z2 - z1);
        // lineVector.Unit(); // Normalize the vector

        // Define a point in the perpendicular plane (e.g., midpoint of the line)
        double midX = (x1 + x2) / 2.0;
        double midY = (y1 + y2) / 2.0;
        double midZ = (z1 + z2) / 2.0;

        // Calculate the normal vector of the perpendicular plane
        // TVector3 normalVector = lineVector.Orthogonal();
        // normalVector.Unit(); // Normalize the vector

        // Calculate the equation of the perpendicular plane: ax + by + cz + d = 0
        // double a = normalVector.X();
        // double b = normalVector.Y();
        // double c = normalVector.Z();
        // double d = -(a * midX + b * midY + c * midZ);

        double a = (x2 - x1);
        double b = (y2 - y1);
        double c = (z2 - z1);
        double d = -(a * midX + b * midY + c * midZ);
        double d1 = -(a * x1 + b * y1 + c * z1);
        double d2 = -(a * x2 + b * y2 + c * z2);
        cout << d1 << " " << d2 << endl;
        double part = histogram->Integral() / 20;
        int div = 0;
        double collect_area = 0;

        if (line < 144)
        {
            div = 200;
            collect_area = 0.1;
        }
        else
        {
            div = 100;
            collect_area = 5;
        }
        for (double ddd = d1; ddd > d2; ddd -= (fabs(d1 - d2)) / div)
        {
            // Print the equation of the perpendicular plane
            // std::cout << "Perpendicular plane equation: " << a << "x + " << b << "y + " << c << "z + " << ddd << " = 0" << std::endl;
            cout << "Counters now " << ddd << ": " << binCount << endl;
            myfile2 << "Counters now " << binCount << " " << line << endl;
            // Count the number of bins in the perpendicular plane
            for (int i = 1; i <= histogram->GetNbinsX() - 100; ++i)
            {
                for (int j = 1; j <= histogram->GetNbinsY() - 100; ++j)
                {
                    for (int k = 1; k <= histogram->GetNbinsZ(); ++k)
                    {
                        if (histogram->GetBinContent(i, j, k) > 0.000000001)
                        {
                            // Get the bin center coordinates
                            double x = histogram->GetXaxis()->GetBinCenter(i);
                            double y = histogram->GetYaxis()->GetBinCenter(j);
                            double z = histogram->GetZaxis()->GetBinCenter(k);
                            // cout << " --------------- " << i << " " << j << " " << k << " coord " << x << " " << y << " " << z << " f(x y z) " << a * x + b * y + c * z + d << " abcd " << a << " " << b << " " << c << " " << d << endl;
                            //  Check if the bin lies in the perpendicular plane
                            if (fabs(a * x + b * y + c * z + ddd) < collect_area)
                            {
                                // cout << "!!!!!!!!!!!!! " << i << " " << j << " " << k << " coord " << x << " " << y << " " << z << " f(x y z) " << a * x + b * y + c * z + d << " abcd " << a << " " << b << " " << c << " " << d << endl;
                                binCount += histogram->GetBinContent(i, j, k);
                                // cout << binCount << " " << i << " " << j << " " << k << " " << histogram->GetBinContent(i, j, k) << endl;
                                if (binCount > 9948 * h_counter)
                                    h_counter++;
                                EdepEmaxMult[h_counter]->SetBinContent(i, j, k, 1); // histogram->GetBinContent(i, j, k)
                                histogram->SetBinContent(i, j, k, 0);

                                // find b_mc and fill
                                vector<double> result = findMatchingLines(i, j, k, h_counter, "data_points_3d_smm_edepemaxmult_b.txt"); // line for edep emax mult b
                                // vector<double> result = findMatchingLines(i, j, k, "data_points_3d_smm_edepemaxmult_npart.txt"); // line for edep emax mult b + npart
                                for (int vv = 0; vv < result.size(); vv++)
                                {
                                    ImpPar[h_counter]->Fill(result[vv]);
                                    // h_npart[h_counter]->Fill(result[vv]);
                                }
                                // cout << i << " " << j << " " << k << " " << result[vv] << " " << endl;
                            }
                        }
                    }
                }
            }
        }
    }

    for (int line = 1; line < 56; line++)
    {
        vector<double> values = read_lines(line, "points_3D_down.txt");
        if (values.size() == 6)
        {
            cout << endl;
            cout << "Line " << line << ": " << values[0] << " " << values[1] << " " << values[2] << endl;
            cout << "Line " << line + 1 << ": " << values[3] << " " << values[4] << " " << values[5] << endl;
        }

        // Define the coordinates of the two points
        double x1 = values[0];
        double y1 = values[1];
        double z1 = values[2];
        double x2 = values[3];
        double y2 = values[4];
        double z2 = values[5];

        // 3.50577 0.630896 7.86
        // 3.60447 0.646994 8.18

        //  Calculate the direction vector of the line
        TVector3 lineVector(x2 - x1, y2 - y1, z2 - z1);
        // lineVector.Unit(); // Normalize the vector

        // Define a point in the perpendicular plane (e.g., midpoint of the line)
        double midX = (x1 + x2) / 2.0;
        double midY = (y1 + y2) / 2.0;
        double midZ = (z1 + z2) / 2.0;

        // Calculate the normal vector of the perpendicular plane
        // TVector3 normalVector = lineVector.Orthogonal();
        // normalVector.Unit(); // Normalize the vector

        // Calculate the equation of the perpendicular plane: ax + by + cz + d = 0
        /*double a = normalVector.X();
        double b = normalVector.Y();
        double c = normalVector.Z();
        double d = -(a * midX + b * midY + c * midZ);*/

        double a = (x2 - x1);
        double b = (y2 - y1);
        double c = (z2 - z1);
        double d = -(a * midX + b * midY + c * midZ);
        double d1 = -(a * x1 + b * y1 + c * z1);
        double d2 = -(a * x2 + b * y2 + c * z2);
        cout << d1 << " " << d2 << endl;
        double part = histogram->Integral() / 20;

        int div = 0;
        double collect_area = 0;

        // div = 50;
        // collect_area = 0.1;

        div = 100;
        collect_area = 5;
        myfile2 << " --------------------2-------------------" << endl;

        for (double ddd = d1; ddd >= d2; ddd -= (fabs(d1 - d2)) / div)
        {
            // Print the equation of the perpendicular plane
            // std::cout << "Perpendicular plane equation: " << a << "x + " << b << "y + " << c << "z + " << ddd << " = 0" << std::endl;
            cout << "Counters now " << binCount << endl;
            myfile2 << "Counters now " << binCount << " " << line << endl;

            // Count the number of bins in the perpendicular plane
            for (int i = 1; i <= histogram->GetNbinsX() - 100; ++i)
            {
                for (int j = 1; j <= histogram->GetNbinsY() - 100; ++j)
                {
                    for (int k = 1; k <= histogram->GetNbinsZ(); ++k)
                    {
                        if (histogram->GetBinContent(i, j, k) > 0.000000001)
                        {
                            // Get the bin center coordinates
                            double x = histogram->GetXaxis()->GetBinCenter(i);
                            double y = histogram->GetYaxis()->GetBinCenter(j);
                            double z = histogram->GetZaxis()->GetBinCenter(k);
                            // cout << " --------------- " << i << " " << j << " " << k << " coord " << x << " " << y << " " << z << " f(x y z) " << a * x + b * y + c * z + d << " abcd " << a << " " << b << " " << c << " " << d << endl;
                            //  Check if the bin lies in the perpendicular plane
                            if (fabs(a * x + b * y + c * z + ddd) < collect_area)
                            {
                                // cout << "!!!!!!!!!!!!! " << i << " " << j << " " << k << " coord " << x << " " << y << " " << z << " f(x y z) " << a * x + b * y + c * z + d << " abcd " << a << " " << b << " " << c << " " << d << endl;
                                binCount += histogram->GetBinContent(i, j, k);

                                if (binCount > 9948 * h_counter)
                                    h_counter++;
                                EdepEmaxMult[h_counter]->SetBinContent(i, j, k, 1);
                                histogram->SetBinContent(i, j, k, 0);
                                // vector<double> result = findMatchingLines(i, j, k, "data_points_3d_smm_edepemaxmult_b.txt");
                                vector<double> result = findMatchingLines(i, j, k, h_counter, "data_points_3d_smm_edepemaxmult_b.txt");
                                for (int vv = 0; vv < result.size(); vv++)
                                    ImpPar[h_counter]->Fill(result[vv]);
                                // h_npart[h_counter]->Fill(result[vv]);
                            }
                        }
                    }
                }
            }
        }
    }
    for (int i = 1; i <= histogram->GetNbinsX() - 100; ++i)
    {
        for (int j = 1; j <= histogram->GetNbinsY() - 100; ++j)
        {
            for (int k = 1; k <= histogram->GetNbinsZ(); ++k)
            {
                if (histogram->GetBinContent(i, j, k) > 0.000000001)
                {
                    vector<double> result = findMatchingLines(i, j, k, h_counter, "data_points_3d_smm_edepemaxmult_b.txt");
                    for (int vv = 0; vv < result.size(); vv++)
                        ImpPar[h_counter]->Fill(result[vv]);
                    // h_npart[h_counter]->Fill(result[vv]);
                }
            }
        }
    }
    // Print the number of bins in the perpendicular plane
    // std::cout << "Bins in the perpendicular plane: " << binCount << std::endl;
    TFile histoFileFull("results_3D_npart_andrey.root", "UPDATE");

    TCanvas *c1 = new TCanvas();
    histogram->Draw();
    histogram->Write();
    cout << "INTEGRAL:  " << histogram->Integral() << endl;

    TCanvas *c2 = new TCanvas();
    for (int i = 0; i < 20; i++)
    {
        if (i != 9)
        {
            EdepEmaxMult[i]->SetMarkerColor(i + 1);
            EdepEmaxMult[i]->Draw("same");
            EdepEmaxMult[i]->Write();
        }
        else
        {
            EdepEmaxMult[i]->SetMarkerColor(5);
            EdepEmaxMult[i]->Draw("same");
            EdepEmaxMult[i]->Write();
        }
    }
    TCanvas *c3 = new TCanvas();
    tmp_for_fit3->Draw();
    tmp_for_fit3->Write();

    TCanvas *c4 = new TCanvas();
    ImpPar[22]->Draw();
    TF1 *fit1 = new TF1("fit1", "gaus");
    for (int i_imp = 0; i_imp < 13; i_imp++)
    {
        cout << "i imp " << i_imp << " integral " << ImpPar[i_imp]->Integral() << endl;
        ImpPar[i_imp]->Fit(fit1, "Q");
        double sigma = fit1->GetParameter(2);
        double mean = fit1->GetParameter(1);
        cout << "mean " << mean << " sigma " << sigma << " s/m " << sigma / mean << endl;
        ImpPar[i_imp]->SetLineColor(i_imp + 1);
        if (ImpPar[i_imp]->Integral() > 0)
        {
            ImpPar[i_imp]->Draw("same");
            ImpPar[i_imp]->Write();
        }
    }

    /*TCanvas *c4 = new TCanvas();
    h_npart[22]->Draw();
    TF1 *fit1 = new TF1("fit1", "gaus");
    for (int i_imp = 0; i_imp < 13; i_imp++)
    {
        cout << "i imp " << i_imp << " integral " << h_npart[i_imp]->Integral() << endl;
        h_npart[i_imp]->Fit(fit1, "Q");
        double sigma = fit1->GetParameter(2);
        double mean = fit1->GetParameter(1);
        cout << "mean fit" << mean << " sigma " << sigma << " s/m " << sigma / mean << " mean " << h_npart[i_imp]->GetMean() << " stddev " << h_npart[i_imp]->GetStdDev() << endl;
        h_npart[i_imp]->SetLineColor(i_imp + 1);
        if (h_npart[i_imp]->Integral() > 0)
        {
            h_npart[i_imp]->Draw("same");
            h_npart[i_imp]->Write();
        }
    }*/

    myfile.close();
    myfile2.close();
}