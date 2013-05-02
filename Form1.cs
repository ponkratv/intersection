using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;

namespace intersection
{
    enum PointClass {LEFT,  RIGHT,  BEYOND,  BEHIND, BETWEEN, ORIGIN, DESTINATION};
    enum CrossingClass {COLLINEAR, PARALLEL, SKEW_NO_CROSS, SKEW_CROSS};
    public partial class Form1 : Form
    {

        Graphics gr;
        public Form1()
        {
            InitializeComponent();
            firstPolygon = new List<Point>();
            secondPolygon = new List<Point>();
            polygons = new List<MyPoint>[2];
            polygons[0] = new List<MyPoint>();
            polygons[1] = new List<MyPoint>();
            gr = panel1.CreateGraphics();
            gr.Clear(Color.White);
        }

        struct MyPoint
        {
            public Point position;
            public bool isIntersection;
            public int indexInNeighbour;
            public bool visited;

            public void setIndexInNeighbour(int val)
            {
                this.indexInNeighbour = val;
            }
        }

        List<Point> firstPolygon, secondPolygon;

        bool selfIntersect(List<Point> poly)
        {
            int n = poly.Count;
            for (int i = 0; i < poly.Count; ++i)
            {
                for (int j = 0; j < poly.Count; ++j)
                {
                    if (i == j) continue;
                    PointF p;
                    if (intersect(poly[i], poly[(i + 1) % n], poly[j], poly[(j + 1) % n], out p))
                    {
                        if (equal(p, poly[i]) || equal(p, poly[(i + 1) % n]))
                        {
                            continue;
                        }
                        return true;
                    }

                }
            }

            return false;
        }
        bool isConvex(List<Point> poly)
        {
            int n = poly.Count;
            int sign = 0;
            for (int i = 0; i < poly.Count; ++i)
            {
                int x1 = poly[(i + 1) % n].X - poly[i].X;
                int y1 = poly[(i + 1) % n].Y - poly[i].Y;

                int x2 = -poly[(i + 1) % n].X + poly[(i + 2) % n].X;
                int y2 = -poly[(i + 1) % n].Y + poly[(i + 2) % n].Y;
                int cross = x1 * y2 - x2 * y1;
                int cur_sign = (cross > 0) ? 1 : -1;
                if (sign == 0)
                {
                    sign = cur_sign;
                }
                else
                {
                    if (sign * cur_sign < 0)
                        return false;
                }
            }
            return true;
        }

        void tryAdd(List<Point> poly, Point p)
        {
            poly.Add(p);
            if (selfIntersect(poly) || !isConvex(poly))
            {
                poly.RemoveAt(poly.Count - 1);
            }
        }

        private void panel1_MouseClick(object sender, MouseEventArgs e)
        {
            if (e.Button == MouseButtons.Left)
            {
                tryAdd(firstPolygon, (new Point(e.X, e.Y)));
                redraw();
            }
            else if (e.Button == MouseButtons.Right)
            {
                tryAdd(secondPolygon, new Point(e.X, e.Y));
                redraw();
            }
        }

        bool inside(double a, double b, double c)
        {
            return b >= Math.Min(a, c) - 1e-5 && b <= Math.Max(a, c) + 1e-5;
        }

        bool intersect(Point a, Point b, Point c, Point d, out PointF p)
        {
            double a1 = a.Y - b.Y;
            double b1 = b.X - a.X;
            double c1 = -(a1 * a.X + b1 * a.Y);

            double a2 = c.Y - d.Y;
            double b2 = d.X - c.X;
            double c2 = -(a2 * c.X + b2 * c.Y);


            PointF res = new Point();
            p = res;
            double dd = a1 * b2 - a2 * b1;
            if (Math.Abs(dd) < 1e-5)
            {
                return false;
            }
            res.X = (float)((b1 * c2 - c1 * b2) / dd);
            res.Y = (float)((a2 * c1 - a1 * c2) / dd);
            if (!inside(a.X, res.X, b.X))
                return false;
            if (!inside(c.X, res.X, d.X))
                return false;
            if (!inside(a.Y, res.Y, b.Y))
                return false;
            if (!inside(c.Y, res.Y, d.Y))
                return false;
            p = res;
            return true;
        }

        bool intersect(PointF p1, PointF p2, double a, double b, double c, out PointF res)
        {
            double pa = p2.Y - p1.Y;
            double pb = p1.X - p2.X;
            double pc = -(p1.X * pa + p1.Y * pb);

            double d = pa * b - a * pb;
            res = new Point();
            if (Math.Abs(d) < 1e-5)
            {
                return false;
            }
            res.X = (float)((pb * c - pc * b) / d);
            res.Y = (float)((a * pc - pa * c) / d);
            if (!inside(p1.X, res.X, p2.X))
                return false;
            if (!inside(p1.Y, res.Y, p2.Y))
                return false;
            return true;
        }

        List<MyPoint>[] polygons;

        double dist(PointF a, PointF b)
        {
            return Math.Sqrt((a.X - b.X) * (a.X - b.X) + (a.Y - b.Y) * (a.Y - b.Y));
        }


        void sort(List<PointF> points, PointF p0)
        {
            double[] angles = new double[points.Count];
            for (int i = 1; i < points.Count; ++i)
            {
                angles[i] = Math.Atan2(points[i].Y - points[0].Y, points[i].X - points[i].X);
            }
            for (int i = 1; i < points.Count; ++i)
            {
                for (int j = i + 1; j < points.Count; ++j)
                {
                    if (angles[i]>angles[j])
                    {
                        PointF p = points[i]; points[i] = points[j]; points[j] = p;
                        double t = angles[i];
                        angles[i] = angles[j];
                        angles[j] = t;
                    }
                }
            }
        }

        bool equal(PointF a, PointF b)
        {
            return Math.Abs(a.X - b.X) < 1e-5 && Math.Abs(a.Y - b.Y) < 1e-5;
        }

        List<PointF> cutLine(List<PointF> polygon, double a, double b, double c)
        {
            List<PointF> ret = new List<PointF>();
            int idx1 = -1;
            int idx2 = -1;
            PointF p1=new Point();
            PointF p2=new Point();
            int n = polygon.Count;
            for (int i = 0; i < polygon.Count; ++i)
            {
                PointF p;
                if (intersect(polygon[i], polygon[(i+1)%n], a, b, c, out p)) {
                    if (idx1 == -1) {idx1 = i;p1=p;}
                    else {idx2 = i;p2=p;}
                }
                if (polygon[i].X * a + polygon[i].Y * b + c < 0)
                {
                    ret.Add(polygon[i]);
                }
            }
            if (idx1 != -1)
            {
                ret.Add(p1);
                ret.Add(p2);
            }

            if (ret.Count > 0)
            {
                sort(ret, ret[0]);
            }
            return ret;
        }

        PointClass classify(PointF p2, PointF p0, PointF p1)
        {
            PointF a = new PointF(p1.X - p0.X, p1.Y - p0.Y);
            PointF b = new PointF(p2.X - p0.X, p2.Y - p0.Y);
            double sa = a.X * b.Y - b.X * a.Y;
            if (sa > 0.0)
                return PointClass.LEFT;
            if (sa < 0.0)
                return PointClass.RIGHT;
            if ((a.X * b.X < 0.0) || (a.Y * b.Y < 0.0))
                return PointClass.BEHIND;
            double lena = Math.Sqrt(a.X * a.X + a.Y * a.Y);
            double lenb = Math.Sqrt(b.X * b.X + b.Y * b.Y);
            if (lena < lenb)
                return PointClass.BEYOND;
            if (p0 == p2)
                return PointClass.ORIGIN;
            if (p1 == p2)
                return PointClass.DESTINATION;
            return PointClass.BETWEEN;
        }

        bool aimsAt (Point a, Point b, Point c, Point d, PointClass aclass, CrossingClass crossType)
        {
            Point va = new Point(b.X - a.X, b.Y - a.Y);
            Point vb = new Point(d.X - c.X, d.Y - c.Y);

            if (crossType != CrossingClass.COLLINEAR)
            {
                if  ( (va.X * vb.Y) >= (vb.X * va.Y) )
                    return (aclass != PointClass.RIGHT);
                else
                    return (aclass != PointClass.LEFT);
            } else {
                return (aclass != PointClass.BEYOND);
            }
        }

        CrossingClass crossingPoint (Point a, Point b, Point c, Point d, out PointF p) {
            int a1 = a.Y - b.Y;
            int b1 = b.X - a.X;
            int c1 = -(a1 * a.X + b1 * a.Y);

            int a2 = c.Y - d.Y;
            int b2 = d.X - c.X;
            int c2 = -(a2 * c.X + b2 * c.Y);

            PointF res = new Point();
            p = res;
            int dd = a1 * b2 - a2 * b1;
            if (dd == 0)
            {
                PointClass aclass = classify(a, c, d);
                if((aclass == PointClass.LEFT) || (aclass == PointClass.RIGHT))
                    return CrossingClass.PARALLEL;
                return CrossingClass.COLLINEAR;
            }
            res.X = (float)((b1 * c2 - c1 * b2) / dd);
            res.Y = (float)((a2 * c1 - a1 * c2) / dd);
            if (!inside(a.X, res.X, b.X))
                return CrossingClass.SKEW_NO_CROSS;
            if (!inside(c.X, res.X, d.X))
                return CrossingClass.SKEW_NO_CROSS;
            if (!inside(a.Y, res.Y, b.Y))
                return CrossingClass.SKEW_NO_CROSS;
            if (!inside(c.Y, res.Y, d.Y))
                return CrossingClass.SKEW_NO_CROSS;
            p = res;
            return CrossingClass.SKEW_CROSS;
        }

        int triangle_area_2(Point a, Point b, Point c)
        {
            return (b.X - a.X) * (c.Y - a.Y) - (b.Y - a.Y) * (c.X - a.X);
        }


        bool pointInConvexPolygon(Point p, List<Point> polygon)
        {
            bool fl1 = false;
            bool fl2 = false;
            for (int i = 0; i < polygon.Count; i++)
            {
                int area = triangle_area_2(polygon[i], polygon[(i + 1) % polygon.Count], p);
                if (area < 0)
                {
                    fl1 = true;
                }
                if (area > 0)
                {
                    fl2 = true;
                }
            }
            if (fl1 && fl2)
                return false;
            return true;
        }

        List<PointF> intersect()
        {
            int n = firstPolygon.Count;
            int m = secondPolygon.Count;

            int inflag = 0;
            int phase = 1;
            int maxIters = 2 * (n + m);

            List<PointF> res = new List<PointF>();
            PointF startPoint = new PointF();
            int i = 0, j = 0;
            for (int k = 0; (k < maxIters) || (phase == 2); k++)
            {
                Point pa = firstPolygon[i];
                Point pb = firstPolygon[(i + n - 1) % n];

                Point qa = secondPolygon[j];
                Point qb = secondPolygon[(j + m - 1) % m];

                PointClass pclass = classify(pb, qa, qb);
                PointClass qclass = classify(qb, pa, pb);

                PointF IntP;
                CrossingClass crossType = crossingPoint(pa, pb, qa, qb, out IntP);

                if (crossType == CrossingClass.SKEW_CROSS)
                {
                    if (phase == 1)
                    {
                        phase = 2;
                        res.Add(IntP);
                        startPoint = IntP;
                    } else if (!equal(res[res.Count - 1], IntP)) {
                        if (equal(startPoint, IntP))
                        {
                            return res;
                        }
                        res.Add(IntP);
                    }
                    if (pclass == PointClass.RIGHT)
                    {
                        inflag = 1;
                    }
                    else if (qclass == PointClass.RIGHT)
                    {
                        inflag = 2;
                    }
                    else
                    {
                        inflag = 0;
                    }
                }
                else if((crossType == CrossingClass.COLLINEAR) && (pclass != PointClass.BEHIND) && (qclass != PointClass.BEHIND))
                {
                    // ? условие проверить
                    inflag = 0;
                }
                bool pAIMSq = aimsAt(pa, pb, qa, qb, pclass, crossType);
                bool qAIMSp = aimsAt(qa, qb, pa, pb, qclass, crossType);
                if (pAIMSq && qAIMSp)
                {
                    if (inflag == 2 || ((inflag == 0) && (pclass == PointClass.LEFT)))
                        i = (n + i - 1) % n;
                    else
                        j = (m + j - 1) % m;
                }
                else if(pAIMSq)
                {
                    i = (n + i - 1) % n;
                    if (inflag == 1 && !equal(res[res.Count - 1], firstPolygon[i]))
                    {
                        res.Add(firstPolygon[i]);
                    }
                }
                else if (qAIMSp)
                {
                    j = (m + j - 1) % m;
                    if (inflag == 2 && !equal(res[res.Count - 1], secondPolygon[j]))
                    {
                        res.Add(secondPolygon[j]);
                    }
                }
                else
                {
                    if ((inflag == 2) || ((inflag == 0) && (pclass == PointClass.LEFT)))
                    {
                        i = (n + i - 1) % n;
                    }
                    else
                    {
                        j = (m + j - 1) % m;
                    }
                }
            }
            res.Clear();
            if (pointInConvexPolygon(firstPolygon[0], secondPolygon))
            {
                for (int k = 0; k < firstPolygon.Count; k++)
                {
                    res.Add(firstPolygon[k]);
                }
            }
            else if(pointInConvexPolygon(secondPolygon[0], firstPolygon))
            {
                for (int k = 0; k < secondPolygon.Count; k++)
                {
                    res.Add(secondPolygon[k]);
                }
            }
            return res;
        }

        bool makeClockWise(List<PointF> p)
        {
            double sq = 0;
            for (int i = 0; i < p.Count; i++)
            {
                PointF p1, p2;
                if (i != 0) p1 = p[i - 1];
                else p1 = p.Last();
                p2 = p[i];
                sq += (p1.X - p2.X) * (p1.Y + p2.Y);
            }
            if (sq > 0)
            {
                p.Reverse();
                return true;
            }
            else
            {
                return false;
            }
        }

        bool makeClockWise(List<Point> p)
        {
            int sq = 0;
            for (int i = 0; i < p.Count; i++)
            {
                Point p1, p2;
                if (i != 0) p1 = p[i - 1];
                else p1 = p.Last();
                p2 = p[i];
                sq += (p1.X - p2.X) * (p1.Y + p2.Y);
            }
            if (sq > 0)
            {
                p.Reverse();
                return true;
            }
            else
            {
                return false;
            }
        }


        void redraw()
        {
            List<Point> firstPolygon = new List<Point>();
            foreach (Point p in this.firstPolygon) firstPolygon.Add(p);
            makeClockWise(firstPolygon);

            List<Point> secondPolygon = new List<Point>();
            foreach (Point p in this.secondPolygon) secondPolygon.Add(p);
            makeClockWise(secondPolygon);


            gr.Clear(Color.White);
            
            if (firstPolygon.Count > 0)
            {
                gr.FillPolygon(new SolidBrush(Color.Red), firstPolygon.ToArray());
            }
            if (secondPolygon.Count > 0)
            {
                gr.FillPolygon(new SolidBrush(Color.Green), secondPolygon.ToArray());
            }
            if (firstPolygon.Count < 3 || secondPolygon.Count < 3)
            {
                return;
            }
            List<PointF> intersection = intersect();

            if (intersection.Count > 0)
            {
                gr.FillPolygon(new SolidBrush(Color.Orange), intersection.ToArray());
            }

        }

        private void button1_Click(object sender, EventArgs e)
        {
            firstPolygon.Clear();
            secondPolygon.Clear();
            redraw();
        }

        private void Form1_Activated(object sender, EventArgs e)
        {
            redraw();
        }
    }
}
