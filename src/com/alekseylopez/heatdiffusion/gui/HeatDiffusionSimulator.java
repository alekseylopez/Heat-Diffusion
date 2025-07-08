package com.alekseylopez.heatdiffusion.gui;

import com.alekseylopez.heatdiffusion.nativebridge.HeatSolver;
import com.alekseylopez.heatdiffusion.nativebridge.HeatSolver.Grid;
import com.alekseylopez.heatdiffusion.nativebridge.HeatSolver.SimulationParams;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.image.BufferedImage;

public class HeatDiffusionSimulator extends JFrame
{
    private Grid grid;
    private SimulationParams params;



    private class HeatMapPanel extends JPanel
    {
        private BufferedImage heatMapImage;
        
        public HeatMapPanel()
        {
            addMouseListener(new MouseAdapter()
            {
                @Override
                public void mouseClicked(MouseEvent e)
                {
                    if(grid != null)
                    {
                        int x = e.getX() * params.nx / getWidth();
                        int y = e.getY() * params.ny / getHeight();
                        
                        if(x >= 0 && x < params.nx && y >= 0 && y < params.ny)
                        {
                            double value = SwingUtilities.isLeftMouseButton(e) ? 10.0 : 0.0;
                            grid.setValue(x, y, value);
                            repaint();
                        }
                    }
                }
            });
        }
        
        @Override
        protected void paintComponent(Graphics g)
        {
            super.paintComponent(g);
            
            if(grid == null)
                return;
            
            int width = getWidth();
            int height = getHeight();
            
            if(heatMapImage == null || heatMapImage.getWidth() != width || heatMapImage.getHeight() != height)
                heatMapImage = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);

            
            double[] data = grid.getData();
            double minVal = Double.MAX_VALUE;
            double maxVal = Double.MIN_VALUE;
            
            // find min and max values
            for(double value : data)
            {
                minVal = Math.min(minVal, value);
                maxVal = Math.max(maxVal, value);
            }
            
            // avoid division by zero
            if(maxVal == minVal)
                maxVal = minVal + 1e-10;
            
            // create heat map
            for(int y = 0; y < height; y++)
            {
                for(int x = 0; x < width; x++)
                {
                    int gridX = x * params.nx / width;
                    int gridY = y * params.ny / height;
                    
                    if(gridX >= params.nx)
                        gridX = params.nx - 1;
                    if(gridY >= params.ny)
                        gridY = params.ny - 1;
                    
                    double value = data[gridX + gridY * params.nx];
                    double normalized = (value - minVal) / (maxVal - minVal);
                    
                    Color color = getHeatColor(normalized);
                    heatMapImage.setRGB(x, y, color.getRGB());
                }
            }
            
            g.drawImage(heatMapImage, 0, 0, null);
            
            // draw grid lines if small enough
            if(params.nx <= 50 && params.ny <= 50)
            {
                g.setColor(Color.GRAY);

                for(int i = 0; i <= params.nx; i++)
                {
                    int x = i * width / params.nx;
                    g.drawLine(x, 0, x, height);
                }

                for(int j = 0; j <= params.ny; j++)
                {
                    int y = j * height / params.ny;
                    g.drawLine(0, y, width, y);
                }
            }
            
            // draw color scale
            drawColorScale(g, width, height, minVal, maxVal);
        }
        
        private Color getHeatColor(double value)
        {
            // clamp value to [0, 1]
            value = Math.max(0.0, Math.min(1.0, value));
            
            // blue to red color map
            if(value < 0.25)
            {
                // black to blue
                double t = value / 0.25;
                return new Color(0, 0, (int) (255 * t));
            } else if(value < 0.5)
            {
                // blue to cyan
                double t = (value - 0.25) / 0.25;
                return new Color(0, (int) (255 * t), 255);
            } else if(value < 0.75)
            {
                // cyan to yellow
                double t = (value - 0.5) / 0.25;
                return new Color((int) (255 * t), 255, (int) (255 * (1 - t)));
            } else
            {
                // yellow to red
                double t = (value - 0.75) / 0.25;
                return new Color(255, (int) (255 * (1 - t)), 0);
            }
        }
        
        private void drawColorScale(Graphics g, int width, int height, double minVal, double maxVal)
        {
            int scaleWidth = 20;
            int scaleHeight = 200;
            int scaleX = width - scaleWidth - 10;
            int scaleY = 10;
            
            // draw color scale
            for(int i = 0; i < scaleHeight; i++)
            {
                double value = 1.0 - (double)i / scaleHeight;
                Color color = getHeatColor(value);
                g.setColor(color);
                g.fillRect(scaleX, scaleY + i, scaleWidth, 1);
            }
            
            // draw scale border
            g.setColor(Color.BLACK);
            g.drawRect(scaleX, scaleY, scaleWidth, scaleHeight);
            
            // draw scale labels
            g.setColor(Color.BLACK);
            g.setFont(new Font("Arial", Font.PLAIN, 10));
            FontMetrics fm = g.getFontMetrics();
            
            String maxLabel = String.format("%.2f", maxVal);
            String minLabel = String.format("%.2f", minVal);
            
            g.drawString(maxLabel, scaleX + scaleWidth + 2, scaleY + fm.getAscent());
            g.drawString(minLabel, scaleX + scaleWidth + 2, scaleY + scaleHeight);
        }
    }

    @Override
    public void dispose()
    {
        if(grid != null)
            grid.destroy();

        super.dispose();
    }
}
