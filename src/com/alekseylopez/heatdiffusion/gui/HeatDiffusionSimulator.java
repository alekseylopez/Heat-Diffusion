package com.alekseylopez.heatdiffusion.gui;

import com.alekseylopez.heatdiffusion.nativebridge.HeatSolver.Grid;
import com.alekseylopez.heatdiffusion.nativebridge.HeatSolver.SimulationParams;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.image.BufferedImage;

public class HeatDiffusionSimulator extends JFrame
{
    private Grid grid;
    private SimulationParams params;
    private Timer animationTimer;
    private HeatMapPanel heatMapPanel;
    private JPanel controlPanel;
    private JButton startButton, pauseButton, resetButton, stepButton;
    private JSlider dtSlider, alphaSlider, speedSlider;
    private JCheckBox implicitCheckbox;
    private JLabel statusLabel, timeLabel;
    private JTextField nxField, nyField, maxStepsField;
    
    private int currentStep = 0;
    private double currentTime = 0.0;
    private boolean isRunning = false;
    private boolean isPaused = false;

    public HeatDiffusionSimulator()
    {
        params = new SimulationParams();

        initializeGUI();
        resetSimulation();
    }

    private void initializeGUI()
    {
        setTitle("Heat Diffusion Simulator");
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setLayout(new BorderLayout());
        
        // create heat map panel
        heatMapPanel = new HeatMapPanel();
        heatMapPanel.setPreferredSize(new Dimension(600, 600));
        add(heatMapPanel, BorderLayout.CENTER);
        
        // create control panel
        createControlPanel();
        add(controlPanel, BorderLayout.EAST);
        
        // create status panel
        createStatusPanel();
        add(createStatusPanel(), BorderLayout.SOUTH);
        
        // animation timer
        animationTimer = new Timer(50, e ->
        {
            if(isRunning && !isPaused && currentStep < params.maxSteps)
            {
                performTimeStep();
                currentStep++;
                currentTime += params.dt;
                updateStatus();
                heatMapPanel.repaint();
                
                if(currentStep >= params.maxSteps)
                    stopSimulation();
            }
        });
        
        pack();
        setLocationRelativeTo(null);
    }

    private void createControlPanel()
    {
        controlPanel = new JPanel();
        controlPanel.setLayout(new BoxLayout(controlPanel, BoxLayout.Y_AXIS));
        controlPanel.setBorder(BorderFactory.createTitledBorder("Controls"));
        controlPanel.setPreferredSize(new Dimension(300, 600));
        
        // simulation controls
        JPanel simControlPanel = new JPanel(new GridLayout(2, 2, 5, 5));
        simControlPanel.setBorder(BorderFactory.createTitledBorder("Simulation"));
        
        startButton = new JButton("Start");
        startButton.addActionListener(e -> startSimulation());
        
        pauseButton = new JButton("Pause");
        pauseButton.addActionListener(e -> pauseSimulation());
        pauseButton.setEnabled(false);
        
        stepButton = new JButton("Step");
        stepButton.addActionListener(e -> performSingleStep());
        
        resetButton = new JButton("Reset");
        resetButton.addActionListener(e -> resetSimulation());
        
        simControlPanel.add(startButton);
        simControlPanel.add(pauseButton);
        simControlPanel.add(stepButton);
        simControlPanel.add(resetButton);
        
        controlPanel.add(simControlPanel);
        
        // parameters panel
        JPanel paramPanel = new JPanel(new GridBagLayout());
        paramPanel.setBorder(BorderFactory.createTitledBorder("Parameters"));
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.insets = new Insets(2, 2, 2, 2);
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.weightx = 1.0;
        
        // grid size
        gbc.gridx = 0; gbc.gridy = 0;
        paramPanel.add(new JLabel("Grid NX:"), gbc);
        gbc.gridx = 1;
        nxField = new JTextField(String.valueOf(params.nx), 6);
        paramPanel.add(nxField, gbc);
        
        gbc.gridx = 0; gbc.gridy = 1;
        paramPanel.add(new JLabel("Grid NY:"), gbc);
        gbc.gridx = 1;
        nyField = new JTextField(String.valueOf(params.ny), 6);
        paramPanel.add(nyField, gbc);
        
        // max steps
        gbc.gridx = 0; gbc.gridy = 2;
        paramPanel.add(new JLabel("Max Steps:"), gbc);
        gbc.gridx = 1;
        maxStepsField = new JTextField(String.valueOf(params.maxSteps), 6);
        paramPanel.add(maxStepsField, gbc);
        
        // time step slider
        gbc.gridx = 0; gbc.gridy = 3;
        paramPanel.add(new JLabel("Time Step:"), gbc);
        gbc.gridx = 1;
        dtSlider = new JSlider(1, 100, (int) (params.dt * 1000));
        dtSlider.addChangeListener(e ->
        {
            params.dt = dtSlider.getValue() / 1000.0;
            updateStatus();
        });
        paramPanel.add(dtSlider, gbc);
        
        // alpha slider
        gbc.gridx = 0; gbc.gridy = 4;
        paramPanel.add(new JLabel("Alpha:"), gbc);
        gbc.gridx = 1;
        alphaSlider = new JSlider(1, 100, (int) (params.alpha * 100));
        alphaSlider.addChangeListener(e ->
        {
            params.alpha = alphaSlider.getValue() / 100.0;
            updateStatus();
        });
        paramPanel.add(alphaSlider, gbc);
        
        // speed slider
        gbc.gridx = 0; gbc.gridy = 5;
        paramPanel.add(new JLabel("Speed:"), gbc);
        gbc.gridx = 1;
        speedSlider = new JSlider(1, 10, 5);
        speedSlider.addChangeListener(e ->
        {
            int delay = 110 - speedSlider.getValue() * 10;
            animationTimer.setDelay(delay);
        });
        paramPanel.add(speedSlider, gbc);
        
        // implicit checkbox
        gbc.gridx = 0; gbc.gridy = 6;
        gbc.gridwidth = 2;
        implicitCheckbox = new JCheckBox("Use Implicit Method", params.useImplicit);
        implicitCheckbox.addActionListener(e -> params.useImplicit = implicitCheckbox.isSelected());
        paramPanel.add(implicitCheckbox, gbc);
        
        controlPanel.add(paramPanel);
        
        // instructions
        JPanel instructionsPanel = new JPanel();
        instructionsPanel.setBorder(BorderFactory.createTitledBorder("Instructions"));
        instructionsPanel.setLayout(new BoxLayout(instructionsPanel, BoxLayout.Y_AXIS));
        
        JTextArea instructions = new JTextArea
        (
            "• Click on the heat map to add heat sources\n" +
            "• Use sliders to adjust parameters\n" +
            "• Start/Pause to control simulation\n" +
            "• Step for single time step\n" +
            "• Reset to restart simulation"
        );
        instructions.setEditable(false);
        instructions.setOpaque(false);
        instructions.setFont(instructions.getFont().deriveFont(12f));
        instructions.setLineWrap(true);
        instructions.setWrapStyleWord(true);
        
        instructionsPanel.add(instructions);
        controlPanel.add(instructionsPanel);
    }

    private JPanel createStatusPanel()
    {
        JPanel statusPanel = new JPanel(new BorderLayout());
        statusPanel.setBorder(BorderFactory.createLoweredBevelBorder());
        
        statusLabel = new JLabel("Ready");
        timeLabel = new JLabel("Time: 0.000");
        
        statusPanel.add(statusLabel, BorderLayout.WEST);
        statusPanel.add(timeLabel, BorderLayout.EAST);
        
        return statusPanel;
    }

    private void startSimulation()
    {
        if(isPaused)
            isPaused = false;
        else
            isRunning = true;
        
        startButton.setEnabled(false);
        pauseButton.setEnabled(true);
        stepButton.setEnabled(false);
        
        animationTimer.start();
        updateStatus();
    }

    private void pauseSimulation()
    {
        isPaused = true;
        animationTimer.stop();
        
        startButton.setEnabled(true);
        pauseButton.setEnabled(false);
        stepButton.setEnabled(true);
        
        updateStatus();
    }

    private void stopSimulation()
    {
        isRunning = false;
        isPaused = false;
        animationTimer.stop();
        
        startButton.setEnabled(true);
        pauseButton.setEnabled(false);
        stepButton.setEnabled(true);
        
        updateStatus();
    }

    private void resetSimulation()
    {
        stopSimulation();
        
        // update parameters from GUI
        try
        {
            params.nx = Integer.parseInt(nxField.getText());
            params.ny = Integer.parseInt(nyField.getText());
            params.maxSteps = Integer.parseInt(maxStepsField.getText());
        } catch(NumberFormatException e)
        {
            JOptionPane.showMessageDialog(this, "Invalid number format in parameters", "Error", JOptionPane.ERROR_MESSAGE);
            return;
        }
        
        // destroy old grid and create new one
        if(grid != null)
            grid.destroy();
        
        grid = new Grid(params.nx, params.ny, params.dx, params.dy, params.alpha);
        
        // set initial heat source at center
        grid.setValue(params.nx / 2, params.ny / 2, 10.0);
        
        currentStep = 0;
        currentTime = 0.0;
        
        heatMapPanel.repaint();
        updateStatus();
    }

    private void performSingleStep()
    {
        if(currentStep < params.maxSteps)
        {
            performTimeStep();
            currentStep++;
            currentTime += params.dt;
            updateStatus();
            heatMapPanel.repaint();
        }
    }

    private void performTimeStep()
    {
        if(grid != null)
        {
            if(params.useImplicit)
                grid.implicitStep(params.dt, params.theta);
            else
                grid.explicitStep(params.dt);
        }
    }

    private void updateStatus()
    {
        if(isRunning && !isPaused)
            statusLabel.setText("Running...");
        else if (isPaused)
            statusLabel.setText("Paused");
        else
            statusLabel.setText("Ready");
        
        timeLabel.setText(String.format("Time: %.3f, Step: %d, dt: %.3f, α: %.2f", currentTime, currentStep, params.dt, params.alpha));
    }

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
