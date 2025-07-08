package com.alekseylopez.heatdiffusion;

import com.alekseylopez.heatdiffusion.gui.HeatDiffusionSimulator;

import javax.swing.*;
import javax.swing.plaf.DimensionUIResource;
import javax.swing.plaf.FontUIResource;

public class Main
{
    public static void main(String[] args)
    {
        try
        {
            UIManager.setLookAndFeel(UIManager.getLookAndFeel());
        } catch(Exception e)
        {
            System.err.println("Could not set system look and feel: " + e.getMessage());
        }
        
        try
        {
            System.loadLibrary("heat_solver");
            System.out.println("Native library loaded successfully");
        } catch(UnsatisfiedLinkError e)
        {
            showLibraryError(e);
            return;
        }
        
        SwingUtilities.invokeLater(() ->
        {
            try
            {
                HeatDiffusionSimulator simulator = new HeatDiffusionSimulator();
                simulator.setVisible(true);
            } catch(Exception e)
            {
                e.printStackTrace();
                JOptionPane.showMessageDialog(null, "Error starting application: " + e.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
            }
        });
    }

    private static void showLibraryError(UnsatisfiedLinkError e) {
        String message = "Failed to load native library 'heat_solver'\n\n" +
                        "Error: " + e.getMessage() + "\n\n" +
                        "Please ensure that:\n" +
                        "• The native library is compiled for your platform\n" +
                        "• heat_solver.dll (Windows) or libheat_solver.so (Linux/Mac) is in your PATH\n" +
                        "• Or add the library location to java.library.path\n\n" +
                        "To compile the native library:\n" +
                        "Windows: gcc -shared -o heat_solver.dll -I\"%JAVA_HOME%\\include\" -I\"%JAVA_HOME%\\include\\win32\" heat_solver.c\n" +
                        "Linux:   gcc -shared -fPIC -o libheat_solver.so -I\"$JAVA_HOME/include\" -I\"$JAVA_HOME/include/linux\" heat_solver.c\n" +
                        "Mac:     gcc -shared -fPIC -o libheat_solver.dylib -I\"$JAVA_HOME/include\" -I\"$JAVA_HOME/include/darwin\" heat_solver.c";
        
        JTextArea textArea = new JTextArea(message);
        textArea.setEditable(false);
        textArea.setFont(new FontUIResource(FontUIResource.MONOSPACED, FontUIResource.PLAIN, 12));
        textArea.setCaretPosition(0);
        
        JScrollPane scrollPane = new JScrollPane(textArea);
        scrollPane.setPreferredSize(new DimensionUIResource(600, 400));
        
        JOptionPane.showMessageDialog(null, scrollPane, 
            "Native Library Error", JOptionPane.ERROR_MESSAGE);
    }
}