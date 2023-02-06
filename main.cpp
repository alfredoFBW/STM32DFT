/*
 * Using the nucleo board to compute FFT and show it on an LCD screen
 */

 /* 
  * Each sample takes 13 us. I've checked it experimentally,
  * therefore the max sampling rate is 1/13 MHz = 76931 KHz 
  * If a sampling rate Fs = 48 KHz is needed then => Ts = 20.83 us
  * therefore there must be a delay between each sample. The delay
  * will be Ts - T_sadc = 20.83 - 13 us \approx 8 us. 
  * With this delay then Fs = 1/(13 + 8 us) = 47619 Hz \approx 48kHz
  */

/* 
 * The library (TextLCD) has some errors due to the fact that some 
 * functions has been removed from MBED 6. What I did was to replace
 * all "wait_ms()" with "ThisThread::sleep_for()". There could be more
 * errors related to the casting, to solve them just make the conversion 
 * from one type to another explicit.
 */


#include "AnalogIn.h"
#include "InterruptIn.h"
#include "PinNames.h"
#include "ThisThread.h"
#include "mbed.h"
#include "mbed_rtx_conf.h"
#include "mbed_thread.h"
#include "mbed_wait_api.h"
#include "stm32f4xx_hal_dma.h"
#include "TextLCD.h"
#include <cmath>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <ctime>

/* LCD pins */
#define D4_LCD D4
#define D5_LCD D5
#define D6_LCD D6 
#define D7_LCD D7
#define RS D8
#define RW D9
#define EN D10

#define SIGNAL   A5
#define DFT_SIZE 256    /* Sweet spot */
#define F_SAMP   48000           
#define PI       3.14159

/* LCD */
#define NUM_UDCS     8
#define NUM_PIX_VERT 16
#define NUM_PIX_HOR  8
#define NUM_COLS_LCD 16
#define NUM_ROWS_LCD 2

/* MODES */
#define NUM_MODES  4
#define MODE_4kHz  0
#define MODE_8kHz  1
#define MODE_16kHz 2
#define MODE_24kHz 3

AnalogIn ain(SIGNAL);
InterruptIn push_button(BUTTON1);

Timer  t;

/* LCD pins */
TextLCD LCD(RS, EN, D4_LCD, D5_LCD, D6_LCD, D7_LCD);

/* LCD scree User Defined Characters (UDC) (udc_0 = udc_bar[0])*/
const char udc_bar[NUM_UDCS][NUM_PIX_HOR] = 
                { {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xFF}, /* 1 Bar */
                  {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xFF, 0xFF}, /* 2 Bars*/
                  {0x00, 0x00, 0x00, 0x00, 0x00, 0xFF, 0xFF, 0xFF},
                  {0x00, 0x00, 0x00, 0x00, 0xFF, 0xFF, 0xFF, 0xFF},
                  {0x00, 0x00, 0x00, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF},
                  {0x00, 0x00, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF},
                  {0x00, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF},
                  {0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF},

                };

/* Since it can be interrupted we need to tell the compiler*/
volatile int mode = MODE_8kHz;


/* Push button Interrupt Service Routine */
void pb_ISR()
{
        mode = (mode + 1) % NUM_MODES;
}

/* Using Radix 4 DIT. For 256 samples it takes 41 ms */
void dft_calc_radix4(float *dft_abs, float *signal_in)
{
        /* Not neccesary the use of vectors. Just update
         * the DFT in each "k" iteration */
        int k, n;
        float Fp_real, Fp_imag;
        float Fpp_real, Fpp_imag;
        float Gp_real, Gp_imag;
        float Gpp_real, Gpp_imag;
        for(k = 0; k < DFT_SIZE/4; k++){
                Fp_real  = Fp_imag  = 0;
                Fpp_real = Fpp_imag = 0;
                Gp_real  = Gp_imag  = 0;
                Gpp_real = Gpp_imag = 0;
                for(n = 0; n < DFT_SIZE/4; n++){
                        Fp_real  += signal_in[4*n]*cosf(2*PI*k*n/(DFT_SIZE/4));
                        Fp_imag  += signal_in[4*n]*sinf(2*PI*k*n/(DFT_SIZE/4));
                        Fpp_real += signal_in[4*n + 2]*cosf(2*PI*k*(4*n + 2)/DFT_SIZE);
                        Fpp_imag += signal_in[4*n + 2]*sinf(2*PI*k*(4*n + 2)/DFT_SIZE);
                        Gp_real  += signal_in[4*n + 1]*cosf(2*PI*k*(4*n + 1)/DFT_SIZE);
                        Gp_imag  += signal_in[4*n + 1]*sinf(2*PI*k*(4*n + 1)/DFT_SIZE);
                        Gpp_real += signal_in[4*n + 3]*cosf(2*PI*k*(4*n + 3)/DFT_SIZE);
                        Gpp_imag += signal_in[4*n + 3]*sinf(2*PI*k*(4*n + 3)/DFT_SIZE);
                }

                dft_abs[k] = sqrtf(powf((Fp_real + Fpp_real) + (Gp_real + Gpp_real), 2) +
                                  powf((Fp_imag + Fpp_imag) + (Gp_imag + Gpp_imag), 2));

                dft_abs[k+DFT_SIZE/4] = sqrtf(powf((Fp_real - Fpp_real) - (Gp_imag - Gpp_imag), 2) +
                                  powf((Fp_imag - Fpp_imag) - (Gpp_real - Gp_real), 2));

                dft_abs[k+2*DFT_SIZE/4] = sqrtf(powf((Fp_real + Fpp_real) - (Gp_real + Gpp_real), 2) +
                                  powf((Fp_imag + Fpp_imag) - (Gp_imag + Gpp_imag), 2));

                dft_abs[k+3*DFT_SIZE/4] = sqrtf(powf((Fp_real - Fpp_real) + (Gp_imag - Gpp_imag), 2) +
                                  powf((Fp_imag - Fpp_imag) + (- Gp_real + Gpp_real), 2));

        }
}

/* Audio freqs, not interested on DC component */
void signal_DC_removal(float *signal)
{
        int k;
        float mean = 0;
        for(k = 0; k < DFT_SIZE; k++)
                mean += signal[k];
        mean /= DFT_SIZE;
        for(k = 0; k < DFT_SIZE; k++)
                signal[k] -= mean;
}

/* Normalize DFT between 0 -> NUM_PIX_VERT */
void dft_normalize_pix(float *dft_abs)
{
        int k; float max = 0.0;
        for( k = 0; k < DFT_SIZE; k++){
                if(dft_abs[k] > max)
                        max = dft_abs[k];
        }
        for(k = 0; k < DFT_SIZE; k++)
                dft_abs[k]  = (float)(NUM_COLS_LCD*dft_abs[k]/max);
}

/* Remember C division (int) is floor() */

/* 
 * How many samples per LCD. Care, if you round to the next integer then 
 * show in LCD that the max_freq for that mode is
 * max_freq = samples_per_lcd_char()*NUM_COLS_LCD*F_SAMP/DFT_SIZE 
 */
int samples_per_lcd_char()
{
        if(mode == MODE_4kHz)
                return ((int)ceilf((4000.0/NUM_COLS_LCD)/(F_SAMP/DFT_SIZE)));
        else if(mode == MODE_8kHz)
                return ((int)ceilf((8000.0/NUM_COLS_LCD)/(F_SAMP/DFT_SIZE)));
        else if(mode == MODE_16kHz)
                return ((int)ceilf((16000.0/NUM_COLS_LCD)/(F_SAMP/DFT_SIZE)));
        else 
                return ((int)floorf((F_SAMP/2/NUM_COLS_LCD)/(F_SAMP/DFT_SIZE)));         
                /* Not aliasing that's why floor is needed */  
}


/* Save UDCs in LCD screen */
void set_udcs()
{
        unsigned char k;
        for( k = 0; k < NUM_UDCS; k++)
                LCD.setUDC(k, (char *)udc_bar[k]);
}

/* 
 * Show the bar related to its value.
 * Value between 0 and 16, col between 1 and 16
 */
void print_bar(int value, int col)
{
        LCD.locate(col, 1);
        if(value <= NUM_PIX_VERT/2){
                LCD.locate(col, 1);
                if(value > 0)
                        LCD.putc(value - 1);
                else
                        LCD.putc(0);
        }
        else{
                LCD.locate(col, 0);                     /* Locate at top */
                LCD.putc(value - NUM_UDCS - 1);         /* Put top BAR */
                LCD.locate(col, 1);                     /* Locate at bot */
                LCD.putc(7);                            /* Full Bar */
        }
}

/* Shows current mode */
void show_mode(int frequency)
{
        LCD.cls();
        LCD.printf("MODE: %dHz", (int)frequency);
        ThisThread::sleep_for(500ms);      
}


void show_dft(float *dft_abs)
{
        LCD.cls();
        int samp_char = samples_per_lcd_char();
        int k, p;
        float height_value;
        dft_normalize_pix(dft_abs);
        for(k = 0; k < NUM_COLS_LCD; k++){
                p = samp_char*k;
                height_value = 0;
                for(; p < (k+1)*samp_char; p++)
                        height_value += dft_abs[p];
                height_value /= samp_char;               /* Mean + Better visualization */
                print_bar((int)height_value, k);
        }

}


int main()
{
        /* Get DFT_SIZE samples. My uC has Floating Point Unit  */


        int k, now_mode;
        float dft_abs[DFT_SIZE];
        float signal_in[DFT_SIZE];

        /* Interrupt */
        push_button.fall(&pb_ISR);

        LCD.printf("Initializing....");
        LCD.setCursor(TextLCD::CurOff_BlkOff);
        ThisThread::sleep_for(1s);         
        set_udcs();

        now_mode = mode - 1;

        while(1){

                /* Take DFT_SIZE samples. It takes DFT_SIZE/Fs seconds */
                for(k = 0; k < DFT_SIZE; k++){       
                        signal_in[k] = ain.read();      /* It takes 13 us */
                        wait_us(8);                     /* It takes 8 us */
                        /* Each sample takes 21 us => Ts = 21 us => Fs = 48000 Hz */
                }

                /* Minimum frequency resolution is F_SAMP/DFT_SIZE */

                /* We do not want DC */
                signal_DC_removal(signal_in);

                /* It takes 41 ms (256 samples) */
                dft_calc_radix4(dft_abs, signal_in);



                /* Check if mode has changed */
                if(now_mode != mode){
                        now_mode = mode;
                        show_mode(samples_per_lcd_char()*NUM_COLS_LCD*F_SAMP/DFT_SIZE);
                }

                show_dft(dft_abs);


                
        }
        return 0;
}


/* Code ending here */
