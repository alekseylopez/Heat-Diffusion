package com.alekseylopez.heatdiffusion.nativebridge;

public class HeatSolver
{
    static
    {
        try
        {
            System.loadLibrary("heat_solver");
        } catch(UnsatisfiedLinkError e)
        {
            System.err.println("Failed to load native library: " + e.getMessage());
            System.err.println("Make sure heat_solver.dll (Windows) or libheat_solver.so (Linux/Mac) is in java.library.path");
        }
    }
    
    // grid management
    public native long createGrid(int nx, int ny, double dx, double dy, double alpha);
    public native void destroyGrid(long gridPtr);
    public native void setGridValue(long gridPtr, int i, int j, double value);
    public native double getGridValue(long gridPtr, int i, int j);
    public native double[] getGridData(long gridPtr);
    
    // time stepping
    public native void explicitStep(long gridPtr, double dt);
    public native void implicitStep(long gridPtr, double dt, double theta);
    
    // grid wrapper class for easier Java usage
    public static class Grid
    {
        private long nativePtr;
        private int nx, ny;
        private double dx, dy, alpha;
        
        public Grid(int nx, int ny, double dx, double dy, double alpha)
        {
            this.nx = nx;
            this.ny = ny;
            this.dx = dx;
            this.dy = dy;
            this.alpha = alpha;
            
            HeatSolver solver = new HeatSolver();
            this.nativePtr = solver.createGrid(nx, ny, dx, dy, alpha);
        }
        
        public void destroy()
        {
            if(nativePtr != 0)
            {
                HeatSolver solver = new HeatSolver();
                solver.destroyGrid(nativePtr);
                nativePtr = 0;
            }
        }
        
        public void setValue(int i, int j, double value)
        {
            HeatSolver solver = new HeatSolver();
            solver.setGridValue(nativePtr, i, j, value);
        }
        
        public double getValue(int i, int j)
        {
            HeatSolver solver = new HeatSolver();
            return solver.getGridValue(nativePtr, i, j);
        }
        
        public double[] getData()
        {
            HeatSolver solver = new HeatSolver();
            return solver.getGridData(nativePtr);
        }
        
        public void explicitStep(double dt)
        {
            HeatSolver solver = new HeatSolver();
            solver.explicitStep(nativePtr, dt);
        }
        
        public void implicitStep(double dt, double theta)
        {
            HeatSolver solver = new HeatSolver();
            solver.implicitStep(nativePtr, dt, theta);
        }
        
        public int getNx()
        {
            return nx;
        }

        public int getNy()
        {
            return ny;
        }
        
        public double getDx()
        {
            return dx;
        }

        public double getDy()
        {
            return dy;
        }

        public double getAlpha()
        {
            return alpha;
        }
    }
    
    // simulation parameters class
    public static class SimulationParams
    {
        public int nx = 100;
        public int ny = 100;
        public double dx = 1.0;
        public double dy = 1.0;
        public double alpha = 0.1;
        public double dt = 0.01;
        public int maxSteps = 1000;
        public boolean useImplicit = true;
        public double theta = 1.0; // 1.0 for backward Euler, 0.5 for Crank-Nicolson
        
        public SimulationParams() {}
        
        public SimulationParams(int nx, int ny, double dx, double dy, double alpha, double dt, int maxSteps, boolean useImplicit, double theta)
        {
            this.nx = nx;
            this.ny = ny;
            this.dx = dx;
            this.dy = dy;
            this.alpha = alpha;
            this.dt = dt;
            this.maxSteps = maxSteps;
            this.useImplicit = useImplicit;
            this.theta = theta;
        }
    }
}