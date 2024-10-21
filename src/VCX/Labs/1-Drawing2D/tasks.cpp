#include <random>
#include <algorithm>
#include <random>

#include <spdlog/spdlog.h>

#include <iostream>

#include "Labs/1-Drawing2D/tasks.h"

using namespace std;

using VCX::Labs::Common::ImageRGB;

namespace VCX::Labs::Drawing2D {
    /******************* 1.Image Dithering *****************/
    void DitheringThreshold(
        ImageRGB &       output,
        ImageRGB const & input) {
        for (std::size_t x = 0; x < input.GetSizeX(); ++x)
            for (std::size_t y = 0; y < input.GetSizeY(); ++y) {
                glm::vec3 color = input.At(x, y);
                output.At(x, y) = {
                    color.r > 0.5 ? 1 : 0,
                    color.g > 0.5 ? 1 : 0,
                    color.b > 0.5 ? 1 : 0,
                };
            }
    }

    void DitheringRandomUniform(
        ImageRGB &       output,
        ImageRGB const & input) {
        // your code here:
        std::mt19937 gen(123);
        uniform_real_distribution<> dist(-0.5,0.5);
        for (std::size_t x = 0; x < input.GetSizeX(); ++x)
            for (std::size_t y = 0; y < input.GetSizeY(); ++y) {
                glm::vec3 color = input.At(x, y);
                double noise = dist(gen);
                output.At(x, y) = {
                    color.r+noise > 0.5 ? 1 : 0,
                    color.g+noise > 0.5 ? 1 : 0,
                    color.b+noise > 0.5 ? 1 : 0,
                };
            }
    }

    void DitheringRandomBlueNoise(
        ImageRGB &       output,
        ImageRGB const & input,
        ImageRGB const & noise) {
            std::mt19937 gen(123);
            uniform_real_distribution<> dist(-1,1);
            ImageRGB tmp(input.GetSizeX(),input.GetSizeY());
            for (std::size_t x = 0; x < input.GetSizeX(); ++x)
            for (std::size_t y = 0; y < input.GetSizeY(); ++y) {
                glm::vec3 color = input.At(x, y);
                glm::vec3 color2 = noise.At(x,y);
                double tval=dist(gen);
                tmp.At(x, y) = {
                    color.r+color2.r-0.5,
                    color.g+color2.g-0.5,
                    color.b+color2.b-0.5,
                };
            }
            DitheringThreshold(output,tmp);
        // your code here:
    }

    void DitheringOrdered(
        ImageRGB &       output,
        ImageRGB const & input) {
            const int array[3][3]={{6,8,4},{1,0,3},{5,2,7}};
            for (std::size_t x = 0; x < input.GetSizeX(); ++x)
            for (std::size_t y = 0; y < input.GetSizeY(); ++y) {
                glm::vec3 color = input.At(x, y);
                for(int i=0;i<=2;++i)for(int j=0;j<=2;++j)
                {
                    int posx=3*x+i,posy=3*y+j;
                    output.At(posx,posy)={
                        color.r*9>array[i][j]?1:0,
                        color.g*9>array[i][j]?1:0,
                        color.b*9>array[i][j]?1:0,
                    };
                }
            }
        // your code here:
    }

    void DitheringErrorDiffuse(
        ImageRGB &       output,
        ImageRGB const & input) {
            ImageRGB tinput=input;
            const int dirx[]={0,1,1,1};
            const int diry[]={1,-1,0,1};
            const double pron[]={7.0/16,3.0/16,5.0/16,1.0/16};
            for (std::size_t x = 0; x < tinput.GetSizeX(); ++x)
            for (std::size_t y = 0; y < tinput.GetSizeY(); ++y) {
                glm::vec3 color = tinput.At(x, y);
                glm::vec3 ans = {
                    color.r > 0.5 ? 1 : 0,
                    color.g > 0.5 ? 1 : 0,
                    color.b > 0.5 ? 1 : 0,
                };
                output.At(x,y)=ans;
                glm::vec3 duff = {
                    color.r-ans.r,
                    color.g-ans.g,
                    color.b-ans.b,
                };
                for(int i=0;i<=3;++i)
                {
                    int nx=x+dirx[i];
                    int ny=y+diry[i];
                    double now=pron[i];
                    if(nx>=0&&nx<(int)tinput.GetSizeX()&&ny>=0&&ny<tinput.GetSizeY())
                    {
                        glm::vec3 tmp=tinput.At(nx,ny);
                        tinput.At(nx,ny)={   
                            tmp.r+now*duff.r,
                            tmp.g+now*duff.g,
                            tmp.b+now*duff.b,
                        };
                    }
                }
            }
        // your code here:
    }

    /******************* 2.Image Filtering *****************/
    void Blur(
        ImageRGB &       output,
        ImageRGB const & input) {
        // your code here:
        const int val[3][3]={{1,1,1},{1,1,1},{1,1,1}};
        for (std::size_t x = 0; x < input.GetSizeX(); ++x)
            for (std::size_t y = 0; y < input.GetSizeY(); ++y) {
                double now=0;
                glm:: vec3 ans={0,0,0};
                int cnt=0;bool flag=0;
                for(int i=-1;i<=1;++i)for(int j=-1;j<=1;++j)
                {
                    int nowx=i+x,nowy=i+y;
                    if(nowx<0||nowx>=input.GetSizeX()||nowy<0||nowy>=input.GetSizeY())
                    {
                        flag=1;
                        break;
                    }
                    glm::vec3 tmp=input.At(nowx,nowy);
                    ans={
                        ans.r+val[i+1][j+1]*tmp.r/9,
                        ans.g+val[i+1][j+1]*tmp.g/9,
                        ans.b+val[i+1][j+1]*tmp.b/9,
                    };
                }
                if(flag==1)output.At(x,y)=input.At(x,y);
                else output.At(x,y)=ans;
            }
    }

    void Edge(
        ImageRGB &       output,
        ImageRGB const & input) {
        // your code here:
        const int val1[3][3]={{-1,0,1},{-2,0,2},{-1,0,1}};
        const int val2[3][3]={{1,2,1},{0,0,0},{-1,-2,-1}};
        for (std::size_t x = 0; x < input.GetSizeX(); ++x)
            for (std::size_t y = 0; y < input.GetSizeY(); ++y) {
                double now=0;
                glm:: vec3 ans1={0,0,0},ans2={0,0,0};
                int cnt=0;bool flag=0;
                for(int i=-1;i<=1;++i)for(int j=-1;j<=1;++j)
                {
                    int nowx=i+x,nowy=i+y;
                    if(nowx<0||nowx>=input.GetSizeX()||nowy<0||nowy>=input.GetSizeY())continue;
                    glm::vec3 tmp=input.At(nowx,nowy);
                    ans1={
                        ans1.r+val1[i+1][j+1]*tmp.r,
                        ans1.g+val1[i+1][j+1]*tmp.g,
                        ans1.b+val1[i+1][j+1]*tmp.b,
                    };
                    ans2={
                        ans2.r+val2[i+1][j+1]*tmp.r,
                        ans2.g+val2[i+1][j+1]*tmp.g,
                        ans2.b+val2[i+1][j+1]*tmp.b,
                    };
                }
                if(flag==1)output.At(x,y)={0,0,0};
                else output.At(x,y)={
                  sqrt(ans1.r*ans1.r+ans2.r*ans2.r), 
                  sqrt(ans1.g*ans1.g+ans2.g*ans2.g),  
                  sqrt(ans1.b*ans1.b+ans2.b*ans2.b)
                };
            }
    }

    /******************* 3. Image Inpainting *****************/
    void Inpainting(
        ImageRGB &         output,
        ImageRGB const &   inputBack,
        ImageRGB const &   inputFront,
        const glm::ivec2 & offset) {
        output             = inputBack;
        std::size_t width  = inputFront.GetSizeX();
        std::size_t height = inputFront.GetSizeY();
        glm::vec3 * g      = new glm::vec3[width * height];
        memset(g, 0, sizeof(glm::vec3) * width * height);
        // set boundary condition
        for (std::size_t y = 0; y < height; ++y) {
            g[y*width]=inputBack.At(offset.x,offset.y+y)-inputFront.At(0,y);
            g[y*width+width-1]=inputBack.At(offset.x+width-1,offset.y+y)-inputFront.At(width-1,y);
            // set boundary for (0, y), your code: g[y * width] = ?
            // set boundary for (width - 1, y), your code: g[y * width + width - 1] = ?
        }
        for (std::size_t x = 0; x < width; ++x) {
            g[x]=inputBack.At(offset.x+x,offset.y)-inputFront.At(x,0);
            g[x+(height-1)*width]=inputBack.At(offset.x+x,offset.y+height-1)-inputFront.At(x,height-1);
            // set boundary for (x, 0), your code: g[x] = ?
            // set boundary for (x, height - 1), your code: g[(height - 1) * width + x] = ?
        }

        // Jacobi iteration, solve Ag = b
        for (int iter = 0; iter < 8000; ++iter) {
            for (std::size_t y = 1; y < height - 1; ++y)
                for (std::size_t x = 1; x < width - 1; ++x) {
                    g[y * width + x] = (g[(y - 1) * width + x] + g[(y + 1) * width + x] + g[y * width + x - 1] + g[y * width + x + 1]);
                    g[y * width + x] = g[y * width + x] * glm::vec3(0.25);
                }
        }

        for (std::size_t y = 0; y < inputFront.GetSizeY(); ++y)
            for (std::size_t x = 0; x < inputFront.GetSizeX(); ++x) {
                glm::vec3 color = g[y * width + x] + inputFront.At(x, y);
                output.At(x + offset.x, y + offset.y) = color;
            }
        delete[] g;
    }

    /******************* 4. Line Drawing *****************/
    void DrawLine(
        ImageRGB &       canvas,
        glm::vec3 const  color,
        glm::ivec2 const p0,
        glm::ivec2 const p1) {
            glm::ivec2 p3=p0,p4=p1;
            if(p3.x>p4.x)swap(p3,p4);
            glm::ivec2 dir=p4-p3;
            glm::ivec2 ndir={-dir.y,dir.x};
            if(abs(dir.y)>abs(dir.x))
            {
                if(dir.y>0)
                {
                    int x=p3.x;
                    int dx=2*(p4.x-p3.x),dy=2*(p4.y-p3.y);
                    int dydx=dy-dx,F=dy/2-dx;
                    for(int y=p3.y;y<=p4.y;++y)
                    {
                        canvas.At(x,y)=color;
                        if(F>0)F-=dx;
                        else x++,F+=dydx;
                    }
                }
                else
                {
                    int x=p3.x;
                    int dx=2*(p4.x-p3.x),dy=2*(p4.y-p3.y);
                    int dydx=dy+dx,F=dy/2+dx;
                    for(int y=p3.y;y>=p4.y;--y)
                    {
                        canvas.At(x,y)=color;
                        if(F<0)F+=dx;
                        else x++,F+=dydx;
                    }
                }
            }
            else
            {
                if(dir.y>0)
                {
                    int y=p3.y;
                    int dx=2*(p4.x-p3.x),dy=2*(p4.y-p3.y);
                    int dydx=dy-dx,F=dy-dx/2;
                    for(int x=p3.x;x<=p4.x;++x)
                    {
                        canvas.At(x,y)=color;
                        if(F<0)F+=dy;
                        else y++,F+=dydx;
                    }
                }
                else
                {
                    int y=p3.y;
                    int dx=2*(p4.x-p3.x),dy=2*(p4.y-p3.y);
                    int dydx=dy+dx,F=dy+dx/2;
                    for(int x=p3.x;x<=p4.x;++x)
                    {
                        canvas.At(x,y)=color;
                        if(F>0)F+=dy;
                        else y--,F+=dydx;
                    }
                }
            }
        // your code here:
    }

    /******************* 5. Triangle Drawing *****************/
    void DrawTriangleFilled(
        ImageRGB &       canvas,
        glm::vec3 const  color,
        glm::ivec2 const p0,
        glm::ivec2 const p1,
        glm::ivec2 const p2) {
            glm:: ivec2 q0=p0,q1=p1,q2=p2;
            if(q0.y<q1.y)swap(q0,q1);
            if(q0.y<q2.y)swap(q0,q2);
            if(q1.y<q2.y)swap(q1,q2);
            glm::ivec2 d1=q1-q0,d2=q2-q0,d3=q2-q1;
            double x1=q0.x,x2=q0.x;
            //std::cerr<<q0.x<<' '<<q0.y<<"\n";
            //std::cerr<<q1.x<<' '<<q1.y<<"\n";
            //std::cerr<<q2.x<<' '<<q2.y<<"\n";
            for(int i=q0.y;i>q1.y;--i)
            {
                //std::cerr<<x1<<' '<<x2<<std::endl;
                for(int j=min(x1,x2);j<=max(x1,x2);++j)canvas.At(j,i)=color;
                if(d1.x!=0)x1+=1.0*d1.x/d1.y;
                if(d2.x!=0)x2+=1.0*d2.x/d2.y;
            }
            for(int i=q1.y;i>=q2.y;--i)
            {
                //std::cerr<<x1<<' '<<x2<<std::endl;
                for(int j=min(x1,x2);j<=max(x1,x2);++j)canvas.At(j,i)=color;
                if(d3.x!=0)x1+=1.0*d3.x/d3.y;
                if(d2.x!=0)x2+=1.0*d2.x/d2.y;
            }
        // your code here:
    }

    /******************* 6. Image Supersampling *****************/
    void Supersample(
        ImageRGB &       output,
        ImageRGB const & input,
        int              rate) {
            //std::cerr<<output.GetSizeX()<<' '<<output.GetSizeY()<<"\n";
            //std::cerr<<input.GetSizeX()<<' '<<input.GetSizeY()<<"\n";
            //std::cerr<<rate<<endl;
            int height=output.GetSizeX(),width=output.GetSizeY();
            vector<glm::vec3>v[height][width];
            double t1=0,t2=0;
            double incx=1.0*input.GetSizeX()/height;
            double incy=1.0*input.GetSizeY()/width;
            for(int i=0;i<output.GetSizeX();++i)
            {
                double t3=t1+incx;t2=0;
                for(int j=0;j<output.GetSizeY();++j)
                {
                    double t4=t2+incy;
                    for(int x=t1;x<int(t3);++x)for(int y=t2;y<int(t4);++y)
                    {
                        glm::vec3 color2=input.At(x,y);
                        color2={
                            color2.r/(rate*rate),
                            color2.g/(rate*rate),
                            color2.b/(rate*rate),
                        };
                        v[i][j].push_back(color2);
                    }
                    t2=t4;
                }
                t1=t3;
            }
            std::mt19937 gen(123);
            for(int i=0;i<output.GetSizeX();++i)
            {
                for(int j=0;j<output.GetSizeY();++j)
                {
                    for(int k=0;k<rate*rate;++k)
                    {
                        int pos=gen()%(int)(v[i][j].size());
                        glm::vec3 color=output.At(i,j);
                        output.At(i,j)={
                            color.r+v[i][j][pos].r,
                            color.g+v[i][j][pos].g,
                            color.b+v[i][j][pos].b,
                        };
                    }
                }
            }

        // your code here:
    }

    /******************* 7. Bezier Curve *****************/
    // Note: Please finish the function [DrawLine] before trying this part.
    glm::vec2 CalculateBezierPoint(
        std::span<glm::vec2> points,
        float const          t) {
        // your code here:
        vector<glm::vec2>point;
        for(auto x:points)point.push_back(x);
        while((int)(point.size())>1)
        {
            vector<glm::vec2>point2;
            for(int j=0;j<(int)(point.size())-1;++j)
            {
                auto now=(1-t)*point[j]+t*point[j+1];
                point2.push_back(now);
            }
            swap(point,point2);
        }
        return point[0];
    }
} // namespace VCX::Labs::Drawing2D