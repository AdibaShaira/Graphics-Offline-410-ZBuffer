#include <bits/stdc++.h>
#include <iostream>
#include <stack>
#include <vector>
#include <cstdio>
#include <fstream>
#include <cmath>
#include <cstring>
#include <string>
#include <sstream>
#include <iomanip>
#include "bitmap_image.hpp"
using namespace std;
struct Colour
{
    double R;
    double G;
    double B;
    void print()
    {
        cout << R << " " << G << " " << B << endl;
    }
};
struct Vector_Axis
{
    double x, y, z;
    void print()
    {
        printf("%.6f %.6f  %.6f", x, y, z);

        //cout<<setprecision(6)<<x<<" "<<y<<" "<<z<<endl;
    }
};
struct Triangle
{
    Vector_Axis point[3];
    Colour col;
    void colourAssign()
    {
        col.R = (rand() % 180);
        col.G = (rand() % 250);
        col.B = (rand() % 230);
    }
    void print()
    {
        for (int i = 0; i < 3; i++)
        {
            //cout<<"point no"<<i<<endl;
            point[i].print();
            cout << endl;
        }

        col.print();
    }
};
struct Homo_Coordinates
{
    double x, y, z, w;
    void normalize()
    {
        if (w != 1)
        {
            x = x / w;
            y = y / w;
            z = z / w;
            w = w / w;
        }
    }
    void print()
    {

        printf("%.1f %.1f %.1f %.1f\n", x, y, z, w);
    }
};

Vector_Axis normalize(Vector_Axis v)
{
    double r;
    r = sqrt((v.x * v.x) + (v.y * v.y) + (v.z * v.z));
    v.x = v.x / r;
    v.y = v.y / r;
    v.z = v.z / r;
    return v;
}
struct Matrix_Structure
{

    double mat[4][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};

    void print()
    {
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                cout << mat[i][j] << " ";
            }
            printf("\n");
        }
    }
};
string eye_pos;
string look_pos;
string up_pos;
Vector_Axis eye_vector;
Vector_Axis look_vector;
Vector_Axis up_vector;
Homo_Coordinates foV;
vector<string> separated_strings[200];
vector<int> push_pop;
stack<Matrix_Structure> matrix_stack;
double Top_Y, Left_X, Bottom_Y, Right_X;
bool left_limit_cross;
void Insert_IdentityMatrix()
{
    Matrix_Structure identity;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            if (i == j)
            {
                //cout<<"dhukse"<<endl;
                identity.mat[i][j] = 1;
            }
            else
            {
                identity.mat[i][j] = 0;
            }
        }
    }
    // cout<<"identity print"<<endl;
    //  identity.print();
    matrix_stack.push(identity);
    Matrix_Structure top = matrix_stack.top();
    //  cout<<"shurute stack"<<endl;
    //  top.print();
}
double Dot_Multiplication(Vector_Axis v1, Vector_Axis v2)
{
    double result;
    result = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    return result;
}
Vector_Axis Normal_Multiplication(Vector_Axis v, double d)
{
    Vector_Axis res;
    res.x = v.x * d;
    res.y = v.y * d;
    res.z = v.z * d;
    return res;
}
Vector_Axis Cross_Multiplication(Vector_Axis v1, Vector_Axis v2)
{
    Vector_Axis res;
    res.x = v1.y * v2.z - v1.z * v2.y;
    res.y = v1.z * v2.x - v1.x * v2.z;
    res.z = v1.x * v2.y - v1.y * v2.x;
    return res;
}
Vector_Axis sum(Vector_Axis a, Vector_Axis b, Vector_Axis c)
{
    Vector_Axis result;
    result.x = a.x + b.x + c.x;
    result.y = a.y + b.y + c.y;
    result.z = a.z + b.z + c.z;
    return result;
}
Vector_Axis Rotate(Vector_Axis x, Vector_Axis a, double angle)
{
    Vector_Axis rotated, t, t1, t2;
    double dot;
    angle = (angle / 180.0) * 3.1416;
    //cout<<"angle radiun "<<angle<<endl;
    dot = Dot_Multiplication(x, a);
    dot = (1 - cos(angle)) * dot;
    t = Normal_Multiplication(a, dot);
    //cout<<"2nd"<<endl;
    //t.print();
    t1 = Normal_Multiplication(x, cos(angle));
    //cout<<"angle  "<<angle<<endl;
    //cout<<"cos   "<<ceil(cos(angle))<<endl;
    //cout<<"1st"<<endl;
    //t1.print();
    t2 = Normal_Multiplication((Cross_Multiplication(a, x)), sin(angle));
    //cout<<"3rd"<<endl;
    //t2.print();
    rotated = sum(t, t1, t2);
    //cout<<"rotated"<<endl;
    rotated.print();
    return rotated;
}
Matrix_Structure Matrix_Multiplication(Matrix_Structure m1, Matrix_Structure m2)
{
    Matrix_Structure mul;
    int row = 4;
    int column = 4;
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < row; j++)
        {
            mul.mat[i][j] = 0;
            for (int k = 0; k < row; k++)
            {
                mul.mat[i][j] += m1.mat[i][k] * m2.mat[k][j];
            }
        }
    }
    // cout<<"matrix mul"<<endl;
    // m1.print();
    // m2.print();
    //  mul.print();
    return mul;
}
void StageOne()
{
    Insert_IdentityMatrix();
    ifstream scene;
    ofstream out;
    scene.open("scene.txt");
    if (!scene.is_open())
    {
        cout << "Error while opening" << endl;
        exit(1);
    }
    string lines;
    int count_line = 0;
    out.open("stageOne.txt");
    while (getline(scene, lines))
    {
        stringstream s(lines);
        string tokens;
        while (s >> tokens)
        {
            separated_strings[count_line].push_back(tokens);
        }
        count_line++;
    }
    //cout<<count_line<<"line no"<<endl;
    for (int i = 0; i <= 3; i++)
    {
        if (i == 0)
        {
            eye_vector.x = stod(separated_strings[i][0].c_str());
            eye_vector.y = stod(separated_strings[i][1].c_str());
            eye_vector.z = stod(separated_strings[i][2].c_str());
        }
        else if (i == 1)
        {
            look_vector.x = stod(separated_strings[i][0].c_str());
            look_vector.y = stod(separated_strings[i][1].c_str());
            look_vector.z = stod(separated_strings[i][2].c_str());
        }

        else if (i == 2)
        {
            up_vector.x = stod(separated_strings[i][0].c_str());
            up_vector.y = stod(separated_strings[i][1].c_str());
            up_vector.z = stod(separated_strings[i][2].c_str());
        }
        else if (i == 3)
        {
            foV.x = stod(separated_strings[i][0].c_str());
            foV.y = stod(separated_strings[i][1].c_str());
            foV.z = stod(separated_strings[i][2].c_str());
            foV.w = stod(separated_strings[i][3].c_str());
        }
    }
    //eye_vector.print();
    //look_vector.print();
    //up_vector.print();
    string command_line;

    for (int i = 3; i < count_line; i++)
    {
        for (int j = 0; j < separated_strings[i].size(); j++)
        {
            command_line = separated_strings[i][j];
            if (command_line == "triangle")
            {
                //ct++;
                Homo_Coordinates tri_first, tri_second, tri_third, transformed_triangle;
                Matrix_Structure top = matrix_stack.top();
                //cout<<"triangle"<<endl;
                //top.print();
                tri_first.x = stod(separated_strings[i + 1][0].c_str());
                tri_first.y = stod(separated_strings[i + 1][1].c_str());
                tri_first.z = stod(separated_strings[i + 1][2].c_str());
                tri_first.w = 1;
                transformed_triangle.x = top.mat[0][0] * tri_first.x + top.mat[0][1] * tri_first.y + top.mat[0][2] * tri_first.z + top.mat[0][3] * tri_first.w;
                transformed_triangle.y = top.mat[1][0] * tri_first.x + top.mat[1][1] * tri_first.y + top.mat[1][2] * tri_first.z + top.mat[1][3] * tri_first.w;
                transformed_triangle.z = top.mat[2][0] * tri_first.x + top.mat[2][1] * tri_first.y + top.mat[2][2] * tri_first.z + top.mat[2][3] * tri_first.w;
                transformed_triangle.w = top.mat[3][0] * tri_first.x + top.mat[3][1] * tri_first.y + top.mat[3][2] * tri_first.z + top.mat[3][3] * tri_first.w;
                if (transformed_triangle.w != 1)
                {
                    transformed_triangle.x = transformed_triangle.x / transformed_triangle.w;
                    transformed_triangle.y = transformed_triangle.y / transformed_triangle.w;
                    transformed_triangle.z = transformed_triangle.z / transformed_triangle.w;
                    transformed_triangle.w = transformed_triangle.w / transformed_triangle.w;
                }
                //transformed_triangle.print();
                out << setprecision(6) << fixed << transformed_triangle.x << " " << transformed_triangle.y << " " << transformed_triangle.z << endl;
                tri_second.x = stod(separated_strings[i + 2][0].c_str());
                tri_second.y = stod(separated_strings[i + 2][1].c_str());
                tri_second.z = stod(separated_strings[i + 2][2].c_str());
                tri_second.w = 1;
                transformed_triangle.x = top.mat[0][0] * tri_second.x + top.mat[0][1] * tri_second.y + top.mat[0][2] * tri_second.z + top.mat[0][3] * tri_second.w;
                transformed_triangle.y = top.mat[1][0] * tri_second.x + top.mat[1][1] * tri_second.y + top.mat[1][2] * tri_second.z + top.mat[1][3] * tri_second.w;
                transformed_triangle.z = top.mat[2][0] * tri_second.x + top.mat[2][1] * tri_second.y + top.mat[2][2] * tri_second.z + top.mat[2][3] * tri_second.w;
                transformed_triangle.w = top.mat[3][0] * tri_second.x + top.mat[3][1] * tri_second.y + top.mat[3][2] * tri_second.z + top.mat[3][3] * tri_second.w;
                if (transformed_triangle.w != 1)
                {
                    transformed_triangle.x = transformed_triangle.x / transformed_triangle.w;
                    transformed_triangle.y = transformed_triangle.y / transformed_triangle.w;
                    transformed_triangle.z = transformed_triangle.z / transformed_triangle.w;
                    transformed_triangle.w = transformed_triangle.w / transformed_triangle.w;
                }
                out << setprecision(6) << fixed << transformed_triangle.x << " " << transformed_triangle.y << " " << transformed_triangle.z << endl;

                tri_third.x = stod(separated_strings[i + 3][0].c_str());
                tri_third.y = stod(separated_strings[i + 3][1].c_str());
                tri_third.z = stod(separated_strings[i + 3][2].c_str());
                tri_third.w = 1;
                transformed_triangle.x = top.mat[0][0] * tri_third.x + top.mat[0][1] * tri_third.y + top.mat[0][2] * tri_third.z + top.mat[0][3] * tri_third.w;
                transformed_triangle.y = top.mat[1][0] * tri_third.x + top.mat[1][1] * tri_third.y + top.mat[1][2] * tri_third.z + top.mat[1][3] * tri_third.w;
                transformed_triangle.z = top.mat[2][0] * tri_third.x + top.mat[2][1] * tri_third.y + top.mat[2][2] * tri_third.z + top.mat[2][3] * tri_third.w;
                transformed_triangle.w = top.mat[3][0] * tri_third.x + top.mat[3][1] * tri_third.y + top.mat[3][2] * tri_third.z + top.mat[3][3] * tri_third.w;
                if (transformed_triangle.w != 1)
                {
                    transformed_triangle.x = transformed_triangle.x / transformed_triangle.w;
                    transformed_triangle.y = transformed_triangle.y / transformed_triangle.w;
                    transformed_triangle.z = transformed_triangle.z / transformed_triangle.w;
                    transformed_triangle.w = transformed_triangle.w / transformed_triangle.w;
                }
                out << setprecision(6) << fixed << transformed_triangle.x << " " << transformed_triangle.y << " " << transformed_triangle.z << endl;
                out << endl;
            }
            else if (command_line == "translate")
            {
                push_pop[push_pop.size() - 1]++;
                double tx, ty, tz;
                // cout<<"dhuksi"<<endl;
                tx = stod(separated_strings[i + 1][0].c_str());
                ty = stod(separated_strings[i + 1][1].c_str());
                tz = stod(separated_strings[i + 1][2].c_str());
                //cout<<"read"<<endl;
                Matrix_Structure mat1;
                //mat1.print();
                for (int k = 0; k < 4; k++)
                {
                    for (int l = 0; l < 4; l++)
                    {
                        // cout<<" loop"<<endl;
                        if (k == 0 && l == 3)
                        {
                            mat1.mat[k][l] = tx;
                        }
                        if (k == 1 && l == 3)
                        {
                            mat1.mat[k][l] = ty;
                        }
                        if (k == 2 && l == 3)
                        {
                            mat1.mat[k][l] = tz;
                        }
                        if (k == l)
                        {
                            mat1.mat[k][l] = 1;
                        }
                    }
                }
                //cout<<"translate print"<<endl;
                //mat1.print();
                Matrix_Structure top = matrix_stack.top();
                Matrix_Structure translated_mat = Matrix_Multiplication(top, mat1);
                matrix_stack.push(translated_mat);
            }
            else if (command_line == "scale")
            {
                push_pop[push_pop.size() - 1]++;
                //cout<<"scale"<<endl;
                double tx, ty, tz;
                tx = stod(separated_strings[i + 1][0].c_str());
                ty = stod(separated_strings[i + 1][1].c_str());
                tz = stod(separated_strings[i + 1][2].c_str());
                //cout<<tx<<ty<<tz<<endl;
                Matrix_Structure mat1;
                for (int k = 0; k < 4; k++)
                {
                    for (int l = 0; l < 4; l++)
                    {
                        if (k == 0 && l == 0)
                        {
                            mat1.mat[k][l] = tx;
                        }
                        if (k == 1 && l == 1)
                        {
                            mat1.mat[k][l] = ty;
                        }
                        if (k == 2 && l == 2)
                        {
                            mat1.mat[k][l] = tz;
                        }
                        if (k == 3 && l == 3)
                        {
                            mat1.mat[k][l] = 1;
                        }
                    }
                }
                //cout<<"after scaling value"<<endl;
                //mat1.print();
                Matrix_Structure top = matrix_stack.top();
                top.print();
                Matrix_Structure scaled_mat = Matrix_Multiplication(top, mat1);
                //cout<<"scaled_mat result"<<endl;
                // scaled_mat.print();
                matrix_stack.push(scaled_mat);
            }
            else if (command_line == "rotate")
            {
                push_pop[push_pop.size() - 1]++;
                double angle, ax, ay, az;
                angle = stod(separated_strings[i + 1][0].c_str());
                ax = stod(separated_strings[i + 1][1].c_str());
                ay = stod(separated_strings[i + 1][2].c_str());
                az = stod(separated_strings[i + 1][3].c_str());
                Vector_Axis a, i, j, k;
                a.x = ax;
                a.y = ay;
                a.z = az;
                i.x = 1;
                i.y = 0;
                i.z = 0;
                j.x = 0;
                j.y = 1;
                j.z = 0;
                k.x = 0;
                k.y = 0;
                k.z = 1;
                a = normalize(a);
                //cout<<"a er val"<<a.x<<" "<<a.y<<" "<<a.z<<endl;
                Vector_Axis c1 = Rotate(i, a, angle);
                //cout<<"c1 er val"<<endl;
                //c1.print();
                Vector_Axis c2 = Rotate(j, a, angle);
                Vector_Axis c3 = Rotate(k, a, angle);
                Matrix_Structure mat1;
                mat1.mat[0][0] = c1.x;
                mat1.mat[0][1] = c2.x;
                mat1.mat[0][2] = c3.x;
                mat1.mat[0][3] = 0;
                mat1.mat[1][0] = c1.y;
                mat1.mat[1][1] = c2.y;
                mat1.mat[1][2] = c3.y;
                mat1.mat[1][3] = 0;
                mat1.mat[2][0] = c1.z;
                mat1.mat[2][1] = c2.z;
                mat1.mat[2][2] = c3.z;
                mat1.mat[2][3] = 0;
                mat1.mat[3][0] = 0;
                mat1.mat[3][1] = 0;
                mat1.mat[3][2] = 0;
                mat1.mat[3][3] = 1;
                // cout<<"rotate er val dhukaisi"<<endl;
                //mat1.print();
                //  mat1.mat[4][4]={{c1.x,c2.x,c3.x,0},{c1.y,c2.y,c3.y,0},{c1.z,c2.z,c3.z,0},{0,0,0,1}};
                Matrix_Structure top = matrix_stack.top();
                //top.print();
                Matrix_Structure rotated_mat = Matrix_Multiplication(top, mat1);
                //cout<<"rotated_mat result"<<endl;
                // rotated_mat.print();
                matrix_stack.push(rotated_mat);
            }
            else if (command_line == "push")
            {
                push_pop.push_back(0);
            }
            else if (command_line == "pop")
            {
                if (matrix_stack.size() == 1)
                {
                    continue;
                }
                int pp = push_pop[push_pop.size() - 1];
                if (pp == 1)
                {
                    push_pop.pop_back();
                    matrix_stack.pop();
                }
                else
                {
                    push_pop.pop_back();
                    while (pp > 0)
                    {
                        matrix_stack.pop();
                        pp--;
                    }
                }
            }
        }
    }
}
void StageTwo()
{
    ifstream scene;
    ofstream out;
    scene.open("stageOne.txt");
    out.open("stageTwo.txt");
    if (!scene.is_open())
    {
        cout << "Error while opening" << endl;
        exit(1);
    }
    string lines;
    int count_line = 0;
    vector<string> stage_two_strings[500];
    while (getline(scene, lines))
    {
        stringstream s(lines);
        string tokens;
        while (s >> tokens)
        {
            stage_two_strings[count_line].push_back(tokens);
        }
        count_line++;
    }
    Vector_Axis l, r, u;
    l.x = look_vector.x - eye_vector.x;
    l.y = look_vector.y - eye_vector.y;
    l.z = look_vector.z - eye_vector.z;
    l = normalize(l);
    r = Cross_Multiplication(l, up_vector);
    r = normalize(r);
    u = Cross_Multiplication(r, l);
    Matrix_Structure Transfor_matrix;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            if (i == j)
            {
                Transfor_matrix.mat[i][j] = 1;
            }
        }
    }
    Transfor_matrix.mat[0][3] = -1.0 * eye_vector.x;
    Transfor_matrix.mat[1][3] = -1.0 * eye_vector.y;
    Transfor_matrix.mat[2][3] = -1.0 * eye_vector.z;
    Matrix_Structure Rotation_matrix;
    Rotation_matrix.mat[0][0] = r.x;
    Rotation_matrix.mat[0][1] = r.y;
    Rotation_matrix.mat[0][2] = r.z;
    Rotation_matrix.mat[1][0] = u.x;
    Rotation_matrix.mat[1][1] = u.y;
    Rotation_matrix.mat[1][2] = u.z;
    Rotation_matrix.mat[2][0] = -1.0 * l.x;
    Rotation_matrix.mat[2][1] = -1.0 * l.y;
    Rotation_matrix.mat[2][2] = -1.0 * l.z;
    Rotation_matrix.mat[3][0] = 0;
    Rotation_matrix.mat[3][1] = 0;
    Rotation_matrix.mat[3][2] = 0;
    Rotation_matrix.mat[3][3] = 1;
    Matrix_Structure V;
    V = Matrix_Multiplication(Rotation_matrix, Transfor_matrix);
    //cout<<"count second stage"<<endl;
    // V.print();

    string command_line;
    for (int i = 0; i < count_line; i++)
    {
        int l = count_line;
        if (i == l - 1)
        {
            break;
        }

        if (stage_two_strings[i].size() == 0)
        {
            i++;
            out << endl;
        }

        Homo_Coordinates tri_first, tri_second, tri_third, transformed_triangle;
        Matrix_Structure top = matrix_stack.top();
        //cout<<"triangle"<<endl;
        //top.print();
        tri_first.x = stod(stage_two_strings[i][0].c_str());
        tri_first.y = stod(stage_two_strings[i][1].c_str());
        tri_first.z = stod(stage_two_strings[i][2].c_str());
        tri_first.w = 1;
        //cout<<"stagetwo input"<<endl;
        //tri_first.print();
        transformed_triangle.x = V.mat[0][0] * tri_first.x + V.mat[0][1] * tri_first.y + V.mat[0][2] * tri_first.z + V.mat[0][3] * tri_first.w;
        transformed_triangle.y = V.mat[1][0] * tri_first.x + V.mat[1][1] * tri_first.y + V.mat[1][2] * tri_first.z + V.mat[1][3] * tri_first.w;
        transformed_triangle.z = V.mat[2][0] * tri_first.x + V.mat[2][1] * tri_first.y + V.mat[2][2] * tri_first.z + V.mat[2][3] * tri_first.w;
        transformed_triangle.w = V.mat[3][0] * tri_first.x + V.mat[3][1] * tri_first.y + V.mat[3][2] * tri_first.z + V.mat[3][3] * tri_first.w;
        if (transformed_triangle.w != 1)
        {
            transformed_triangle.x = transformed_triangle.x / transformed_triangle.w;
            transformed_triangle.y = transformed_triangle.y / transformed_triangle.w;
            transformed_triangle.z = transformed_triangle.z / transformed_triangle.w;
            transformed_triangle.w = transformed_triangle.w / transformed_triangle.w;
        }
        //transformed_triangle.print();
        out << setprecision(6) << fixed << transformed_triangle.x << " " << transformed_triangle.y << " " << transformed_triangle.z << endl;

        //cout<<"for loop count line "<<i<<endl;
    }

    //cout<<"done"<<endl;
}
double Top_Clipping_Y(Triangle t)
{
    double num1, num2, num3;
    double max_y;
    num1 = t.point[0].y;
    num2 = t.point[1].y;
    num3 = t.point[2].y;
    if (num1 >= num2 && num1 >= num3)
    {
        max_y = num1;
    }
    else if (num2 >= num1 && num2 >= num3)
    {
        max_y = num2;
    }
    else
    {
        max_y = num3;
    }
    if (max_y > Top_Y)
    {
        max_y = Top_Y;
    }
    return max_y;
}
double Bottom_Clipping_Y(Triangle t)
{
    double num1, num2, num3;
    double min_y;
    num1 = t.point[0].y;
    num2 = t.point[1].y;
    num3 = t.point[2].y;
    if (num1 <= num2 && num1 <= num3)
    {
        min_y = num1;
    }
    else if (num2 <= num1 && num2 <= num3)
    {
        min_y = num2;
    }
    else
    {
        min_y = num3;
    }
    if (min_y < Bottom_Y)
    {
        min_y = Bottom_Y;
    }
    return min_y;
}
double Left_Clipping_X(Triangle t)
{
    double num1, num2, num3;
    double min_x;
    num1 = t.point[0].x;
    num2 = t.point[1].x;
    num3 = t.point[2].x;
    if (num1 <= num2 && num1 <= num3)
    {
        min_x = num1;
    }
    else if (num2 <= num1 && num2 <= num3)
    {
        min_x = num2;
    }
    else
    {
        min_x = num3;
    }
    //cout<<"min_x....."<<min_x<<endl;
    if (min_x > Right_X)
    {
        left_limit_cross = false;
    }
    if (min_x < Left_X)
    {
        min_x = Left_X;
    }
    return min_x;
}
double Right_Clipping_X(Triangle t)
{
    double num1, num2, num3;

    double max_y;
    num1 = t.point[0].x;
    num2 = t.point[1].x;
    num3 = t.point[2].x;
    if (num1 >= num2 && num1 >= num3)
    {
        max_y = num1;
    }
    else if (num2 >= num1 && num2 >= num3)
    {
        max_y = num2;
    }
    else
    {
        max_y = num3;
    }
    if (max_y > Right_X)
    {
        max_y = Right_X;
    }
    return max_y;
}
void StageThree()
{
    double fovy = foV.x;
    fovy = (fovy / 180.0) * 3.1416;
    double aspectRatio = foV.y;
    double near = foV.z;
    double far = foV.w;
    double fovx = fovy * aspectRatio;
    double t = near * tan(fovy / 2);
    double r = near * tan(fovx / 2);
    Matrix_Structure transform_matrix;
    transform_matrix.mat[0][0] = near / r;
    transform_matrix.mat[0][1] = 0;
    transform_matrix.mat[0][2] = 0;
    transform_matrix.mat[0][3] = 0;
    transform_matrix.mat[1][0] = 0;
    transform_matrix.mat[1][1] = near / t;
    transform_matrix.mat[1][2] = 0;
    transform_matrix.mat[1][3] = 0;
    transform_matrix.mat[2][0] = 0;
    transform_matrix.mat[2][1] = 0;
    transform_matrix.mat[2][2] = (-1.0 * (far + near)) / (far - near);
    transform_matrix.mat[2][3] = (-2.0 * far * near) / (far - near);
    transform_matrix.mat[3][0] = 0;
    transform_matrix.mat[3][1] = 0;
    transform_matrix.mat[3][2] = -1.0;
    transform_matrix.mat[3][3] = 0;
    //cout<<"anishaa.............."<<endl;
    //  transform_matrix.print();
    ifstream scene;
    ofstream out;
    scene.open("stageTwo.txt");
    out.open("stageThree.txt");
    //cout<<"StageThree..........................."<<endl;
    if (!scene.is_open())
    {
        cout << "Error while opening" << endl;
        exit(1);
    }
    string lines;
    int count_line = 0;
    vector<string> stage_two_strings[500];
    while (getline(scene, lines))
    {
        stringstream s(lines);
        string tokens;
        while (s >> tokens)
        {
            stage_two_strings[count_line].push_back(tokens);
        }
        count_line++;
    }
    for (int i = 0; i < count_line; i++)
    {

        int l = count_line;
        if (i == l)
        {
            break;
        }

        if (stage_two_strings[i].size() == 0)
        {
            i++;
            out << endl;
        }

        Homo_Coordinates tri_first, tri_second, tri_third, transformed_triangle;
        //Matrix_Structure top = matrix_stack.top();
        //cout<<"triangle"<<endl;
        //top.print();
        tri_first.x = stod(stage_two_strings[i][0].c_str());
        tri_first.y = stod(stage_two_strings[i][1].c_str());
        tri_first.z = stod(stage_two_strings[i][2].c_str());
        tri_first.w = 1;
        // cout<<"stagetwo input"<<endl;
        // tri_first.print();
        transformed_triangle.x = transform_matrix.mat[0][0] * tri_first.x + transform_matrix.mat[0][1] * tri_first.y + transform_matrix.mat[0][2] * tri_first.z + transform_matrix.mat[0][3] * tri_first.w;
        transformed_triangle.y = transform_matrix.mat[1][0] * tri_first.x + transform_matrix.mat[1][1] * tri_first.y + transform_matrix.mat[1][2] * tri_first.z + transform_matrix.mat[1][3] * tri_first.w;
        transformed_triangle.z = transform_matrix.mat[2][0] * tri_first.x + transform_matrix.mat[2][1] * tri_first.y + transform_matrix.mat[2][2] * tri_first.z + transform_matrix.mat[2][3] * tri_first.w;
        transformed_triangle.w = transform_matrix.mat[3][0] * tri_first.x + transform_matrix.mat[3][1] * tri_first.y + transform_matrix.mat[3][2] * tri_first.z + transform_matrix.mat[3][3] * tri_first.w;
        if (transformed_triangle.w != 1)
        {
            transformed_triangle.x = transformed_triangle.x / transformed_triangle.w;
            transformed_triangle.y = transformed_triangle.y / transformed_triangle.w;
            transformed_triangle.z = transformed_triangle.z / transformed_triangle.w;
            transformed_triangle.w = transformed_triangle.w / transformed_triangle.w;
        }
        //transformed_triangle.print();
        out << setprecision(5) << fixed << transformed_triangle.x << " " << transformed_triangle.y << " " << transformed_triangle.z << endl;

        //cout<<"for loop count line "<<i<<endl;
    }

    //cout<<"ber hoise"<<endl;
}
void StageFour()
{
    cout << "hi" << endl;
    int screen_width, screen_height;
    double left_limit_x;
    double bottom_limit_y;
    double front_limit, rear_limit;
    ifstream scene;
    ifstream config;
    scene.open("config.txt");
    vector<string> stage_config_strings[20];
    string lines;
    int count_line = 0;
    while (getline(scene, lines))
    {
        // cout<<"dhukoo"<<endl;
        stringstream s(lines);
        string tokens;
        while (s >> tokens)
        {
            cout << tokens << endl;
            stage_config_strings[count_line].push_back(tokens);
        }
        count_line++;
    }
    screen_width = stod(stage_config_strings[0][0].c_str());
    screen_height = stod(stage_config_strings[0][1].c_str());
    left_limit_x = stod(stage_config_strings[1][0].c_str());
    bottom_limit_y = stod(stage_config_strings[2][0].c_str());
    front_limit = stod(stage_config_strings[3][0].c_str());
    rear_limit = stod(stage_config_strings[3][1].c_str());
    //cout<<screen_width<<" "<<screen_height<<" "<<left_limit_x<<" "<<bottom_limit_y<<" "<<front_limit<<" "<<rear_limit<<" "<<endl;
    ifstream scene1;
    scene1.open("stageThree.txt");
    lines = " ";
    count_line = 0;
    vector<string> stage_three_strings[500];
    while (getline(scene1, lines))
    {
        //cout<<"dhukoo"<<endl;
        stringstream s(lines);
        string tokens;
        while (s >> tokens)
        {
            // cout<<tokens<<endl;
            stage_three_strings[count_line].push_back(tokens);
        }
        count_line++;
    }

    int triangle_number;
    triangle_number = count_line / 4;
    Triangle tri_obj[triangle_number + 5];
    int c = 0;
    for (int i = 0; i < count_line; i = i + 4)
    {
        if (stage_three_strings[i].size() == 0)
        {
            continue;
        }
        for (int j = 0; j < 3; j++)
        {
            tri_obj[c].point[j].x = stod(stage_three_strings[i + j][0].c_str());
            tri_obj[c].point[j].y = stod(stage_three_strings[i + j][1].c_str());
            tri_obj[c].point[j].z = stod(stage_three_strings[i + j][2].c_str());
        }
        c++;
    }
    for (int i = 0; i <= triangle_number; i++)
    {
        //cout<<"kire"<<endl;
        tri_obj[i].colourAssign();
        //tri_obj[i].print();
    }

    double z_buffer[screen_width][screen_height];
    //vector<vector<int>> frame_buffer;
    for (int i = 0; i < screen_width; i++)
    {
        for (int j = 0; j < screen_height; j++)
        {
            z_buffer[i][j] = rear_limit;
        }
    }
    double dx, dy, top_limit_y, right_limit_x;
    right_limit_x = -1.0 * left_limit_x;
    top_limit_y = -1.0 * bottom_limit_y;
    dx = (right_limit_x - left_limit_x) / screen_width;
    dy = (top_limit_y - bottom_limit_y) / screen_height;
    Top_Y = top_limit_y - (dy / 2);
    Left_X = left_limit_x + (dx / 2);
    Bottom_Y = bottom_limit_y + (dy / 2);
    Right_X = right_limit_x - (dx / 2);
    //cout<<"Top Y"<<Top_Y<<endl;
    //cout<<"Left_X"<<Left_X<<endl;
    //cout<<"Right_X"<<Right_X<<endl;
    //cout<<"Bottom_Y"<<Bottom_Y<<endl;
    bitmap_image image(screen_width, screen_height);
    for (int i = 0; i < screen_width; i++)
    {
        for (int j = 0; j < screen_height; j++)
        {
            image.set_pixel(i, j, 0, 0, 0);
        }
    }

    for (int i = 0; i <= triangle_number; i++)
    {
        left_limit_cross = true;
        double topy = Top_Clipping_Y(tri_obj[i]);
        double bottomy = Bottom_Clipping_Y(tri_obj[i]);
        double leftx = Left_Clipping_X(tri_obj[i]);
        double rightx = Right_Clipping_X(tri_obj[i]);
        if (!left_limit_cross)
        {
            continue;
        }
        cout << "Triangle no " << i << endl;
        cout << topy << "----row-----" << bottomy << endl;
        // cout<<"left_limit..."<<Left_X<<"right_limit....."<<Right_X<<endl;
        // cout<<"i no triangle..........."<<i<<endl;
        // cout<<leftx<<"----coloumn----"<<rightx<<endl;
        int row_no_t = round((Top_Y - topy) / dy);
        int lo = round((leftx - Left_X) / dx);
        int ro = round((rightx + Right_X) / dx);
        // cout<<".....leftrow...."<<lo<<"......rightrow....."<<ro<<endl;
        if (row_no_t < 0)
        {
            row_no_t = -1 * row_no_t;
        }
        int row_no_b = round((Bottom_Y + bottomy) / dy);
        if (row_no_b < 0)
        {
            row_no_b = -1 * row_no_b;
        }
        //cout<<"row no top..."<<row_no_t<<"row no bottom"<<row_no_b<<endl;
        int t = min(row_no_t, row_no_b);
        int b = max(row_no_b, row_no_t);
        //cout<<"ki shomosha..............."<<t<<"....bottom....."<<b<<endl;
        for (int j = t; j < b; j++)
        {
            double xa, ys, xb, xc;
            ys = Top_Y - j * dy;
            //cout<<"ys..."<<ys<<endl;
            double x1, x2, x3, y1, y2, y3;
            x1 = tri_obj[i].point[0].x;
            x2 = tri_obj[i].point[1].x;
            x3 = tri_obj[i].point[2].x;
            y1 = tri_obj[i].point[0].y;
            y2 = tri_obj[i].point[1].y;
            y3 = tri_obj[i].point[2].y;
            //cout<<"points.."<<x1<<" "<<x2<<" "<<x3<<" "<<y1<<" "<<y2<<" "<<y3<<endl;
            xa = x1 + ((ys - y1) * (x2 - x1)) / (y2 - y1);
            xb = x1 + ((ys - y1) * (x3 - x1)) / (y3 - y1);
            xc = x2 + ((ys - y2) * (x2 - x3)) / (y2 - y3);
            double za, zb, zc, z1, z2, z3;
            z1 = tri_obj[i].point[0].z;
            z2 = tri_obj[i].point[1].z;
            z3 = tri_obj[i].point[2].z;
            za = z1 + ((ys - y1) * (z2 - z1)) / (y2 - y1);
            zb = z1 + ((ys - y1) * (z3 - z1)) / (y3 - y1);
            zc = z2 + ((ys - y2) * (z2 - z3)) / (y2 - y3);
            double xa_temp, xb_temp, xc_temp, za_temp, zb_temp;
            bool a, b, c;
            a = true;
            b = true;
            c = true;
            if (xa > 1000000 || isnan(xa))
            {

                xa_temp = xb;
                xb_temp = xc;
                za_temp = zb;
                zb_temp = zc;
            }
            else if (xb > 1000000 || isnan(xb))
            {

                xa_temp = xa;
                xb_temp = xc;
                za_temp = za;
                zb_temp = zc;
            }
            else if (xc > 1000000 || isnan(xc))
            {

                xa_temp = xa;
                xb_temp = xb;
                za_temp = za;
                zb_temp = zb;
            }
            else
            {

                if (xa == xb || xb == xc || xc == xa)
                {
                    if (xa == xb)
                    {
                        xa_temp = xa;
                        xb_temp = xc;
                        za_temp = zb;
                        zb_temp = zc;
                    }
                    else if (xb == xc)
                    {
                        xa_temp = xa;
                        xb_temp = xb;
                        za_temp = za;
                        zb_temp = zc;
                    }
                    else
                    {
                        xa_temp = xa;
                        xb_temp = xb;
                        za_temp = za;
                        zb_temp = zb;
                    }
                }
                else
                {
                    if (xa < leftx || xa > rightx)
                    {

                        xa_temp = xb;
                        xb_temp = xc;
                        za_temp = zb;
                        zb_temp = zc;
                    }
                    else if (xb < leftx || xb > rightx)
                    {

                        xa_temp = xa;
                        xb_temp = xc;
                        za_temp = za;
                        zb_temp = zc;
                    }
                    else if (xc < leftx || xc > rightx)
                    {

                        xa_temp = xa;
                        xb_temp = xb;
                        za_temp = za;
                        zb_temp = zb;
                    }
                }
            }
            //cout<<"xa..."<<xa_temp<<"  xb..."<<xb_temp<<endl;
            //cout<<"left x..."<<Left_X<<" right x..."<<Right_X<<endl;
            //cout<<"dx..."<<dx<<endl;
            if (xa_temp < xb_temp)
            {
                if (xa_temp < Left_X)
                {
                    xa_temp = Left_X;
                }
                if (xb_temp > Right_X)
                {
                    xb_temp = Right_X;
                }
            }
            else if (xa_temp > xb_temp)
            {
                if (xa_temp > Right_X)
                {
                    xa_temp = Right_X;
                }
                if (xb_temp < Left_X)
                {
                    xb_temp = Left_X;
                }
            }
            int col_no_l = round((xa_temp - Left_X) / dx);
            int col_no_r = round((xb_temp + Right_X) / dx);
            /*if(j==230){
             cout<<"j val...."<<j<<" col l..."<<col_no_l<<" col r..."<<col_no_r<<endl;
            }*/
            int col_left = min(col_no_l, col_no_r);
            int col_right = max(col_no_l, col_no_r);
            for (int x = col_left; x <= col_right; x++)
            {
                double xp, zp;
                xp = Left_X + x * dx;
                double change;
                zp = za_temp + ((xp - xa_temp) * (zb_temp - za_temp)) / (xb_temp - xa_temp);
                if (zp < z_buffer[j][x] && zp >= front_limit && zp <= rear_limit)
                {
                    z_buffer[j][x] = zp;
                    image.set_pixel(x, j, tri_obj[i].col.R, tri_obj[i].col.G, tri_obj[i].col.B);
                }
            }
        }
        image.save_image("check.bmp");
        //image.clear();
        ofstream out;
        out.open("zbuffer.txt");
        for (int i = 0; i < screen_width; i++)
        {
            for (int j = 0; j < screen_height; j++)
            {
                if (z_buffer[i][j] == rear_limit)
                {


                }
              else{
                    out << " " << z_buffer[i][j];
                }
            }
            out << endl;
        }
        free(z_buffer);
    }
}
void huhu()
{
    bitmap_image image(500, 300);

    for (int i = 0; i < 200; i++)
    {
        for (int j = 0; j < 100; j++)
        {
            image.set_pixel(i, j, 255, 255, 0);
        }
        for (int j = 150; j < 200; j++)
        {
            image.set_pixel(i, j, 0, 0, 255);
        }
    }

    for (int i = 400; i < 450; i++)
    {
        for (int j = 50; j < 150; j++)
        {
            image.set_pixel(i, j, 255, 0, 0);
        }
        for (int j = 200; j < 300; j++)
        {
            image.set_pixel(i, j, 0, 255, 255);
        }
    }

    image.save_image("test.bmp");
    ;
}
int main()
{
    cout << "hi" << endl;
    StageOne();
    cout << "done" << endl;
    StageTwo();
    StageThree();
    StageFour();
    //huhu();
}
