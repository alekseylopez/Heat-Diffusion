package com.alekseylopez.heatdiffusion;

import com.alekseylopez.heatdiffusion.gui.HeatDiffusionSimulator;

import javax.swing.*;

public class Main
{
    public static void main(String[] args)
    {
        SwingUtilities.invokeLater(() ->
        {
            try
            {
                UIManager.setLookAndFeel(UIManager.getLookAndFeel());
            } catch(Exception e)
            {
                e.printStackTrace();
            }
            
            new HeatDiffusionSimulator().setVisible(true);
        });
    }
}