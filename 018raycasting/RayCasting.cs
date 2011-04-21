﻿using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.Windows.Forms;
using OpenTK;
using Rendering;

namespace Rendering
{
  /// <summary>
  /// Test scene for ray-based rendering.
  /// </summary>
  public class RayScene : IRayScene
  {
    public IIntersectable Intersectable
    {
      get;
      set;
    }

    public double[] BackgroundColor
    {
      get;
      set;
    }

    public ICamera Camera
    {
      get;
      set;
    }

    public ICollection<ILightSource> Sources
    {
      get;
      set;
    }

    /// <summary>
    /// Creates default ray-rendering scene.
    /// </summary>
    public RayScene ()
    {
      // scene:
      CSGInnerNode root = new CSGInnerNode( SetOperation.Union );
      Intersectable = root;
      // sphere 1:
      Sphere s = new Sphere();
      s.SetAttribute( PropertyName.COLOR, new double[] { 1.0, 0.6, 0.0 } );
      root.InsertChild( s, Matrix4d.Identity );
      // sphere 2:
      s = new Sphere();
      s.SetAttribute( PropertyName.COLOR, new double[] { 0.2, 0.9, 0.5 } );
      root.InsertChild( s, Matrix4d.CreateTranslation( -2.2, 0.0, 0.0 ) );
      // sphere 3:
      s = new Sphere();
      s.SetAttribute( PropertyName.COLOR, new double[] { 0.1, 0.3, 1.0 } );
      root.InsertChild( s, Matrix4d.CreateTranslation( -4.4, 0.0, 0.0 ) );
      // sphere 4:
      s = new Sphere();
      s.SetAttribute( PropertyName.COLOR, new double[] { 1.0, 0.2, 0.2 } );
      root.InsertChild( s, Matrix4d.CreateTranslation( 2.2, 0.0, 0.0 ) );
      // sphere 5:
      s = new Sphere();
      s.SetAttribute( PropertyName.COLOR, new double[] { 0.1, 0.4, 0.0 } );
      s.SetAttribute( PropertyName.TEXTURE, new CheckerTexture( 80.0, 40.0, new double[] { 1.0, 0.8, 0.2 } ) );
      root.InsertChild( s, Matrix4d.CreateTranslation( 4.4, 0.0, 0.0 ) );

      // background color:
      BackgroundColor = new double[] { 0.0, 0.1, 0.2 };

      // camera:
      Camera = new StaticCamera( new Vector3d( 0.0, 0.0, -10.0 ),
                                 new Vector3d( 0.0, 0.0, 1.0 ), 60.0 );

      // light sources:
      Sources = new LinkedList<ILightSource>();
      Sources.Add( new PointLightSource( new Vector3d( -10.0, 8.0, 3.0 ), 1.0 ) );
    }
  }

  /// <summary>
  /// Correctly implemented set operations..
  /// </summary>
  public class CSGInnerNode : InnerNode
  {
    public CSGInnerNode ( SetOperation op )
      : base( op )
    { }

    /// <summary>
    /// Computes the complete intersection of the given ray with the object. 
    /// </summary>
    /// <param name="p0">Ray origin.</param>
    /// <param name="p1">Ray direction vector.</param>
    /// <returns>Sorted list of intersection records.</returns>
    public override LinkedList<Intersection> Intersect ( Vector3d p0, Vector3d p1 )
    {
      // !!!{{ TODO: put correct set operation implementation here

      LinkedList<Intersection> result = new LinkedList<Intersection>();
      foreach ( ISceneNode n in children )
      {
        Vector3d origin = Vector3d.TransformPosition( p0, n.FromParent );
        Vector3d dir = Vector3d.TransformVector( p1, n.FromParent );
        // ray in local child's coords: [ origin, dir ]

        LinkedList<Intersection> partial = n.Intersect( origin, dir );
        if ( partial == null || partial.Count == 0 ) continue;

        Intersection i = null;
        foreach ( Intersection inter in partial )
          if ( inter.T > 0.0 )
          {
            i = inter;
            break;
          }

        if ( i != null )
          if ( result.First == null )
            result.AddFirst( i );
          else
            if ( i.T < result.First.Value.T )
              result.AddFirst( i );
      }
      return result;

      // !!!}}
    }
  }

}

namespace _018raycasting
{
  public partial class Form1 : Form
  {
    private void redraw ()
    {
      Cursor.Current = Cursors.WaitCursor;

      int width = panel1.Width;
      int height = panel1.Height;
      outputImage = new Bitmap( width, height, System.Drawing.Imaging.PixelFormat.Format24bppRgb );

      SimpleImageSynthesizer sis = new SimpleImageSynthesizer();
      sis.Width = width;
      sis.Height = height;

      // default constructor of the RayScene .. custom scene
      RayScene scene = new RayScene();
      IImageFunction imf = new RayCasting( scene );
      imf.Width = width;
      imf.Height = height;
      sis.ImageFunction = imf;

      Stopwatch sw = new Stopwatch();
      sw.Start();

      sis.RenderRectangle( outputImage, 0, 0, width, height );

      sw.Stop();
      labelElapsed.Text = String.Format( "Elapsed: {0:f}s", 1.0e-3 * sw.ElapsedMilliseconds );

      pictureBox1.Image = outputImage;

      Cursor.Current = Cursors.Default;
    }

  }
}
