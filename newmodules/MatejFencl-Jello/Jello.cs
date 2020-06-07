using MathSupport;
using OpenTK;
using Rendering;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO.Compression;
using System.Xml.Linq;
using Utilities;

namespace MatejFencl_Jello
{
	public class JelloManager
	{
		private static double floorLevel = -1;

		public static void SetFloor(double f)
		{
			floorLevel = f;
		}

		public static double GetFloor ()
		{
			return floorLevel;
		}

		private static Dictionary<int,Jello>  jello = new Dictionary<int, Jello>();

		public static Jello RegisterJello (int id,int XS,int YS,int ZS,double S,double XO, double YO, double ZO,double J=1.0)
		{
			Jello j = new Jello(XS,YS,ZS,S,XO,YO,ZO,J);
			jello[id] = j;
			return j;
		}

		public static Jello GetJello (int id)
		{
			return jello[id];
		}
	}

	public class JelloMesh : DefaultSceneNode, ISolid, ITimeDependent
	{

		public void GetBoundingBox (out Vector3d min, out Vector3d max)
		{
			min = Vector3d.Zero;
			max = Vector3d.Zero;
		}

		public int Frame = 0;

		public override LinkedList<Intersection> Intersect (Vector3d p0, Vector3d p1)
		{

			if (Frame < 0 || Frame > jello.Frames)
				return null;

			List<Intersection> result = new List<Intersection>();


			var isect1 = surf_top.Intersect(p0, p1);
			if(isect1 != null)
				result.AddRange(isect1);
			var isect2=surf_bottom.Intersect(p0, p1);
			if (isect2 != null)
				result.AddRange(isect2);
			var isect3=surf_left.Intersect(p0, p1);
			if (isect3 != null)
				result.AddRange(isect3);
			var isect4=surf_right.Intersect(p0, p1);
			if (isect4 != null)
				result.AddRange(isect4);
			var isect5=surf_front.Intersect(p0, p1);
			if (isect5 != null)
				result.AddRange(isect5);
			var isect6=surf_back.Intersect(p0, p1);
			if (isect6 != null)
				result.AddRange(isect6);
			result.Sort();

			foreach (Intersection isect in result)
			{
				isect.Solid = this;
			}

			return new LinkedList<Intersection>(result);
		}


		public override void CompleteIntersection (Intersection inter)
		{
			if (surf_top!=null)
				surf_top.CompleteIntersection(inter);
		}


		public int getSerial ()
		{
			return 0;
		}

		public object Clone ()
		{
			return new JelloMesh(this);
		}

		BezierSurface surf_top,surf_bottom,surf_left, surf_right,surf_front,surf_back;

		void update_in_time()
		{
			Frame = jello.getFrame(time);
			
			surf_top = null;
			surf_bottom = null;
			surf_left = null;
			surf_right = null;
			surf_front = null;
			surf_back = null;

			System.GC.Collect();

			if (Frame < 0 || Frame > jello.Frames)
				return;

			surf_top = new BezierSurface((jello.ZS-1)/3, (jello.XS - 1) / 3, jello.vertices_top[Frame]);
			surf_bottom = new BezierSurface((jello.ZS - 1) / 3, (jello.XS - 1) / 3, jello.vertices_bottom[Frame]);
			surf_left = new BezierSurface((jello.YS - 1) / 3, (jello.ZS - 1)/3, jello.vertices_left[Frame]);
			surf_right = new BezierSurface((jello.YS - 1) / 3, (jello.ZS - 1) / 3, jello.vertices_right[Frame]);
			surf_front = new BezierSurface((jello.YS - 1) / 3, (jello.XS - 1) / 3, jello.vertices_front[Frame]);
			surf_back = new BezierSurface((jello.YS - 1) / 3, (jello.XS - 1) / 3, jello.vertices_back[Frame]);

			surf_top.PreciseNormals = true;
			surf_bottom.PreciseNormals = true;
			surf_left.PreciseNormals = true;
			surf_right.PreciseNormals = true;
			surf_front.PreciseNormals = true;
			surf_back.PreciseNormals = true;
		}

		Jello jello;

		public double Start { get { return 0; } set { } }
		public double End { get { return double.PositiveInfinity; } set { } }

		double time;
		public double Time { get { return time; } set { time = value; update_in_time(); } }

		public JelloMesh (int id)
		{
			jello = JelloManager.GetJello(id);
		}

		public JelloMesh (JelloMesh jm)
		{
			jello = jm.jello;
		}
	}

	[Serializable]
	public class Jello
	{
		public class VolumePoint
		{
			public double fac;
			public double x, y, z;
			public double xv, yv, zv;
		}

		public int XS,YS,ZS;
		public double J;

		public double XO,YO,ZO;
		public double S;

		public Jello (int xs, int ys, int zs, double s, double xo,double yo,double zo, double j = 1.0)
		{
			XS = xs*3+1;
			YS = ys*3+1;
			ZS = zs*3+1;
			XO = xo;
			YO = yo;
			ZO = zo;
			J = j;
			S = s;
		}

		public void AddForce(double XO, double YO, double ZO, double XD, double YD, double ZD, double radius)
		{
			default_forces.Add(new Tuple<Vector3, Vector3, double>(new Vector3((float)XO, (float)YO, (float)ZO), new Vector3((float)XD, (float)YD, (float)ZD), radius));
		}

		public List<double[]> vertices_top = new List<double[]>();
		public List<double[]> vertices_bottom = new List<double[]>();
		public List<double[]> vertices_left = new List<double[]>();
		public List<double[]> vertices_right = new List<double[]>();
		public List<double[]> vertices_front = new List<double[]>();
		public List<double[]> vertices_back = new List<double[]>();

		public List<Tuple<Vector3,Vector3,double>> default_forces = new List<Tuple<Vector3, Vector3, double>>();

		public int Frames = 0;
		public double FPS = 0;
		public double StartTime = 0;
		public double EndTime = 0;

		public int getFrame (double time)
		{
			return (int)Math.Round((time - StartTime) / (EndTime - StartTime) * Frames);
		}

		void React (VolumePoint p1, VolumePoint p2, double req_dist, double fac)
		{
			double force = (Math.Sqrt((p1.x-p2.x)*(p1.x-p2.x) +  (p1.y-p2.y)*(p1.y-p2.y) +(p1.z-p2.z)*(p1.z-p2.z)) - req_dist)/ fac;

			force /= J;

			Vector3 n;
			n.X = (float)(p2.x - p1.x);
			n.Y = (float)(p2.y - p1.y);
			n.Z = (float)(p2.z - p1.z);
			n.Normalize();

			p1.xv += n.X * force;
			p1.yv += n.Y * force;
			p1.zv += n.Z * force;
			p2.xv -= n.X * force;
			p2.yv -= n.Y * force;
			p2.zv -= n.Z * force;
		}

		public void Simulate (double StartTime, double EndTime, double FPS)
		{
			this.StartTime = StartTime;
			this.EndTime = EndTime;
			this.FPS = FPS;
			this.Frames = (int)Math.Round((EndTime - StartTime) * FPS);

			List<Tuple<VolumePoint,VolumePoint,double>> links = new List<Tuple<VolumePoint, VolumePoint, double>>();


			for (int f = 0; f < Frames; f++)
			{
				vertices_top.Add(new double[XS*ZS*3]);
				vertices_bottom.Add(new double[XS*ZS*3]);
				vertices_left.Add(new double[ZS*YS*3]);
				vertices_right.Add(new double[ZS*YS*3]);
				vertices_front.Add(new double[XS*YS*3]);
				vertices_back.Add(new double[XS*YS*3]);
			}


			VolumePoint[,,] volumePoints = new VolumePoint[XS,YS,ZS];


			for (int x = 0; x < XS; x++)
			{
				for (int y = 0; y < YS; y++)
				{
					for (int z = 0; z < ZS; z++)
					{
						VolumePoint vp = new VolumePoint();

						vp.x = XO+((double)(x - (XS - 1) / 2.0) / 3.0) * S;
						vp.y = YO+((double)(y - (YS - 1) / 2.0) / 3.0) * S;
						vp.z = ZO+((double)(z - (ZS - 1) / 2.0) / 3.0) * S;

						vp.xv = 0;
						vp.yv = 0;
						vp.zv = 0;

						foreach(var force in default_forces)
						{
							double d = Math.Sqrt((vp.x - force.Item1.X)*(vp.x - force.Item1.X) + (vp.y - force.Item1.Y)*(vp.y - force.Item1.Y) + (vp.z - force.Item1.Z)*(vp.z - force.Item1.Z));
							double f = Math.Max(0,force.Item3 - d) / force.Item3 / 100.0;
							vp.xv += force.Item2.X * f;
							vp.yv += force.Item2.Y * f;
							vp.zv += force.Item2.Z * f;
						}

						vp.fac = 1.0;
						volumePoints[x, y, z] = vp;
					}
				}
			}

			for (int x = 0; x < XS; x++)
			{
				for (int y = 0; y < YS; y++)
				{
					for (int z = 0; z < ZS; z++)
					{
						for (int xa = -2; xa <= 2; xa++)
						{
							for (int ya = -2; ya <= 2; ya++)
							{
								for (int za = -2; za <= 2; za++)
								{
									if (xa != 0 && ya != 0 && za != 0 && x + xa > 0 && y + ya > 0 && z + za > 0 && x + xa < XS && y + ya < YS && z + za < ZS)
									{
										VolumePoint vp1 = volumePoints[x, y, z];
										VolumePoint vp2 = volumePoints[x + xa, y + ya, z + za];
										vp1.fac += 1;
										double xd = vp2.x - vp1.x;
										double yd = vp2.y - vp1.y;
										double zd = vp2.z - vp1.z;
										double dist = Math.Sqrt(xd* xd + yd*yd + zd*zd);
										links.Add(new Tuple<VolumePoint, VolumePoint, double>(vp1,vp2,  dist));
									}
								}
							}
						}
					}
				}
			}

			for (int f = 0; f < Frames; f++)
			{
				for (int i = 0; i < vertices_top[f].Length; i++)
					vertices_top[f][i] = 0.0;
				for (int i = 0; i < vertices_bottom[f].Length; i++)
					vertices_bottom[f][i] = 0.0;
				for (int i = 0; i < vertices_left[f].Length; i++)
					vertices_left[f][i] = 0.0;
				for (int i = 0; i < vertices_right[f].Length; i++)
					vertices_right[f][i] = 0.0;
				for (int i = 0; i < vertices_front[f].Length; i++)
					vertices_front[f][i] = 0.0;
				for (int i = 0; i < vertices_back[f].Length; i++)
					vertices_back[f][i] = 0.0;

				for (int step = 0; step < 10; step++)
				{
					for (int x = 0; x < XS; x++)
					{
						for (int y = 0; y < YS; y++)
						{
							for (int z = 0; z < ZS; z++)
							{
								VolumePoint p = volumePoints[x,y,z];
								p.yv -= 0.00004;
							}
						}
					}

					foreach (var link in links)
					{
						React(link.Item1, link.Item2, link.Item3, link.Item1.fac);
					}

					for (int x = 0; x < XS; x++)
					{
						for (int y = 0; y < YS; y++)
						{
							for (int z = 0; z < ZS; z++)
							{
								VolumePoint p = volumePoints[x,y,z];
								p.x += p.xv;
								p.y += p.yv;
								p.z += p.zv;

								p.xv *= 0.998;
								p.yv *= 0.998;
								p.zv *= 0.998;

								if (p.y < JelloManager.GetFloor())
								{
									p.yv += JelloManager.GetFloor() - p.y;
									p.xv /= 2.0;
									p.zv /= 2.0;
								}
							}
						}
					}



					for (int x = 0; x < XS; x++)
					{
						for (int z = 0; z < ZS; z++)
						{
							VolumePoint vp = volumePoints[x, 0,z];
							VolumePoint vp1 = volumePoints[x,YS-1,z];


							int id = x+z*XS;
							vertices_top[f][id * 3] += vp.x;
							vertices_top[f][id * 3 + 1] += vp.y;
							vertices_top[f][id * 3 + 2] += vp.z;

							int id1 = (XS-1-x)+z*XS;
							vertices_bottom[f][id1 * 3] += vp1.x;
							vertices_bottom[f][id1 * 3 + 1] += vp1.y;
							vertices_bottom[f][id1 * 3 + 2] += vp1.z;
						}
					}

					for (int z = 0; z < ZS; z++)
					{
						for (int y = 0; y < YS; y++)
						{
							VolumePoint vp = volumePoints[0, y,z];
							VolumePoint vp1 = volumePoints[XS-1,y,z];


							int id = z+y*ZS;
							vertices_left[f][id * 3] += vp.x;
							vertices_left[f][id * 3 + 1] += vp.y;
							vertices_left[f][id * 3 + 2] += vp.z;

							int id1 = (ZS-1-z)+y*ZS;
							vertices_right[f][id1 * 3] += vp1.x;
							vertices_right[f][id1 * 3 + 1] += vp1.y;
							vertices_right[f][id1 * 3 + 2] += vp1.z;
						}
					}

					for (int x = 0; x < XS; x++)
					{
						for (int y = 0; y < YS; y++)
						{
							VolumePoint vp = volumePoints[x, y,0];
							VolumePoint vp1 = volumePoints[x,y,ZS-1];

							int id = (XS-1-x)+y*XS;
							vertices_front[f][id * 3] += vp.x;
							vertices_front[f][id * 3 + 1] += vp.y;
							vertices_front[f][id * 3 + 2] += vp.z;

							int id1 = x+y*XS;
							vertices_back[f][id1 * 3] += vp1.x;
							vertices_back[f][id1 * 3 + 1] += vp1.y;
							vertices_back[f][id1 * 3 + 2] += vp1.z;
						}
					}

				}

				for (int i = 0; i < vertices_top[f].Length; i++)
					vertices_top[f][i] /= 10.0;
				for (int i = 0; i < vertices_bottom[f].Length; i++)
					vertices_bottom[f][i] /= 10.0;
				for (int i = 0; i < vertices_left[f].Length; i++)
					vertices_left[f][i] /= 10.0;
				for (int i = 0; i < vertices_right[f].Length; i++)
					vertices_right[f][i] /= 10.0;
				for (int i = 0; i < vertices_front[f].Length; i++)
					vertices_front[f][i] /= 10.0;
				for (int i = 0; i < vertices_back[f].Length; i++)
					vertices_back[f][i] /= 10.0;


			}
		}


	}
}
