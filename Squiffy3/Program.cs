using System;
using System.IO;

namespace Squiffy3
{
    namespace SlantScanBasic
    {
        class Program
        {
            //const int hRawFrame = 128;
            static void Main(string[] args)
            {

                Console.WriteLine("Hello World!");

                //lets read in a file and the "gain"
                string folder;
                string fn = "";
                //string fnGain = folder + "slantGain002.tif";
                int ppS = 87;
                int stepScale = 1019;
                // var imgInp = ezTiff.ReadTiffMulti(folder + fn);


                folder = @"C:\Users\Xela\Downloads\";
                fn = @"Result of Substack (500-600)rev.raw";
                var fnGain = folder + @"Gain Result of Substack (500-600)rev.raw";
                //folder = @"C:\Users\spencer\OneDrive\code\Alex\";
                //fn = folder + @"Result of Substack (500-600)rev.raw";
                //fnG = folder + @"Gain Result of Substack (500-600)rev.raw";
                ppS = 95;
                stepScale = 1087;

                BinaryReader br = new BinaryReader(File.OpenRead(folder + fn));
                var imgInpFrames = new Single[100][];
                for (int z = 0; z < 100; z++)
                {

                    imgInpFrames[z] = new float[128 * 1024];

                    for (int n = 0; n < imgInpFrames[z].Length; n++)
                    {
                        imgInpFrames[z][n] = (float)br.ReadByte();

                    }
                }
                br.Close();
                br.Dispose();

                //var imgInpFrames = ezFile.EviFileStuff.ReadRawChunks(
                //    Path.Join(folder, fn),
                //    tc: TypeCode.Byte,
                //    Width: 1024,
                //    nChunks: 100,
                //    linesPerChunk: 128,
                //    chunkSkip: 0,
                //    hdrBytes: 0);
                //    var imgInp = (width: 1024, height: 128, pix: imgInpFrames);
                var imgInpW = 1024;
                const int imgInpH = 128;
                var imgpix = imgInpFrames;
                //    read data in 




                var dimInpw = imgInpW;
                var dimIpnh = imgInpH;
                var w = dimInpw;
                var h = dimIpnh;
                var dimOutW = imgInpW * 4;
                const int yFramePad4 = imgInpH * 4;
                var dimOutH = 2*yFramePad4 + (int)((imgpix.Length * 4L * 1100L) / 1000);
                var dimOut2W = dimOutW / 2;
                var dimOut2H = dimOutH / 2;
                //unpack to img0
                var img0 = new float[dimOutW * dimOutH];
                //weight of img0 pixels
                var imgw = new float[dimOutW * dimOutH];
                //unsheared output
                var img1 = new float[dimOutW * dimOutH];
                //unsheared output weight
                var img1w = new float[dimOutW * dimOutH];
                //Post filtered output
                var img2f = new float[dimOutW * dimOutH];
                //Post filtered output weight
                var img2w = new float[dimOutW * dimOutH];
                //binned image
                var imgb = new float[dimOut2W * dimOut2H];
                //binned output weight
                var imgbw = new float[dimOut2W * dimOut2H];
                //output
                var imgOut = new float[dimOut2W * dimOut2H];

                int w4 = dimOutW;
                int h4 = dimOutH;
                int len4 = h4 * w4;
                int w2 = dimOut2W;
                int h2 = dimOut2H;
                bool UnSlant = true;
                if (UnSlant)
                {

                    ///////////////////////////////////////////////////////////////////
                    //Forward Project and Unpack
                    //
                    for (int nz = 0; nz < imgpix.Length; nz++)
                    {
                        if (nz % 1000 == 0) { Console.Write("*"); } // progress
                        if (nz == 99)
                        {
                            Console.Write("*");
                        } // progress
                        for (int ny = 0; ny < h; ny++)
                        {
                            for (int nx = 0; nx < w; nx++)
                            {
                                //Add blank columns 
                                //- change the source index by dnx which is the number of tile gaps 
                                int dnx = (nx / 129) + 1;//includes the tile gap 
                                if (nx % 129 == 0)
                                {
                                    continue;
                                }
                                //unbin
                                int nx_0 = nx * 4;
                                int ny_0 = ny * 4;
                                //int nz_0 = nz * 4;
                                // Shear X for tilt
                                int nx0 = nx_0 - ((ny_0 * ppS) / 1000);
                                //Motion in y per frame
                                int ny0 = ny_0 - ((nz * 4 * stepScale) / 1000);
                                ny0 += h4 - (yFramePad4*2);

                                int nPix = ny0 * w4 + nx0;
                                //nPix += 1 + w4;
                                if (nPix < 0) { continue; }
                                if (nPix >= len4) { continue; }
                                //if (nPix + 3 * w4 + 3 >= len4) { continue; }

                                int niPix = ny * w + nx - dnx; // ghost pixels
                                float vPix = imgpix[nz][niPix];// * imgGain[niPix];

                                //ALEX Post Blur !
                                //ALEX? weight by gain
                                img0[nPix] += vPix;
                                imgw[nPix] += 1;

                            }
                        }
                    }


                    //Ok now we have unpacked
                    BinaryWriter bw = new BinaryWriter(File.Create(@"slant1.raw"));
                    for (int n = 0; n < img0.Length; n++)
                    {
                        bw.Write(img0[n]);
                    }
                    bw.Close();
                    bw.Dispose();

                    //for (int ny = 0; ny < h4; ny++)
                    //{
                    //    for (int nx = 0; nx < w4; nx++)
                    //    {
                    //        //invert imgw
                    //        int nPix = ny * w4 + nx;
                    //        float t = Math.Max(1, imgw[nPix]);
                    //        imgw[nPix] = 1.0f / t;
                    //    }
                    //}
                    //for (int ny = 0; ny < h4; ny++)
                    //{
                    //    for (int nx = 0; nx < w4; nx++)
                    //    {
                    //        //ALEX FIX
                    //        //apply weight
                    //        int nPix = ny * w4 + nx;
                    //        img0[nPix] *= imgw[nPix];
                    //    }
                    //}

                    ///////////////////////////////////////////////////////////
                    //
                    // Final Unshear
                    //-out side of itterative loop?
                    //
                    //int dh4 = (((w4 * ppS) / 1000));
                    for (int ny = 0; ny < h4; ny++)
                    {
                        for (int nx = 0; nx < w4; nx++)
                        {
                            //Un-Shear columns
                            int nPix1 = ny * w4 + nx;
                            //Whole quaterpixel shear
                            //?dither
                            int nPix0 = (ny - ((nx * ppS) / 1000)) * w4 + nx;
                            if (nPix0 < 0) { continue; }
                            if (nPix0 >= len4) { continue; }
                            img1[nPix1] = img0[nPix0];
                            img1w[nPix1] = imgw[nPix0];
                        }
                    }
                    //PostBlur
                    //ALEX ? Separable !
                    for (int ny = 0; ny < h4; ny++)
                    {
                        for (int nx = 0; nx < w4; nx++)
                        {
                            int nPix4 = ny * w4 + nx;
                            if (nPix4 < 0) { continue; }
                            if (nPix4 + 3 * w4 + 3 >= len4) { continue; }

                            float vPix = img1[nPix4];
                            float wpix = img1w[nPix4];
                            img2f[nPix4] += 2 * vPix;
                            img2w[nPix4] += 2 * wpix;
                            img2f[nPix4 + 1] += 3 * vPix;
                            img2w[nPix4 + 1] += 3 * wpix;
                            img2f[nPix4 + 2] += 3 * vPix;
                            img2w[nPix4 + 2] += 3 * wpix;
                            img2f[nPix4 + 3] += 2 * vPix;
                            img2w[nPix4 + 3] += 2 * wpix;

                            img2f[nPix4 + w4] += 3 * vPix;
                            img2w[nPix4 + w4] += 3 * wpix;
                            img2f[nPix4 + w4 + 1] += 4 * vPix;
                            img2w[nPix4 + w4 + 1] += 4 * wpix;
                            img2f[nPix4 + w4 + 2] += 4 * vPix;
                            img2w[nPix4 + w4 + 2] += 4 * wpix;
                            img2f[nPix4 + w4 + 3] += 3 * vPix;
                            img2w[nPix4 + w4 + 3] += 3 * wpix;

                            img2f[nPix4 + w4 + w4] += 3 * vPix;
                            img2w[nPix4 + w4 + w4] += 3 * wpix;
                            img2f[nPix4 + w4 + w4 + 1] += 4 * vPix;
                            img2w[nPix4 + w4 + w4 + 1] += 4 * wpix;
                            img2f[nPix4 + w4 + w4 + 2] += 4 * vPix;
                            img2w[nPix4 + w4 + w4 + 2] += 4 * wpix;
                            img2f[nPix4 + w4 + w4 + 3] += 3 * vPix;
                            img2w[nPix4 + w4 + w4 + 3] += 3 * wpix;

                            img2f[nPix4 + w4 + w4 + w4] += 2 * vPix;
                            img2w[nPix4 + w4 + w4 + w4] += 2 * wpix;
                            img2f[nPix4 + w4 + w4 + w4 + 1] += 3 * vPix;
                            img2w[nPix4 + w4 + w4 + w4 + 1] += 3 * wpix;
                            img2f[nPix4 + w4 + w4 + w4 + 2] += 3 * vPix;
                            img2w[nPix4 + w4 + w4 + w4 + 2] += 3 * wpix;
                            img2f[nPix4 + w4 + w4 + w4 + 3] += 2 * vPix;
                            img2w[nPix4 + w4 + w4 + w4 + 3] += 2 * wpix;
                        }
                    }
                    img1 = img2f;
                    img1w = img2w;
                    //weighted bin
                    //ALEX DO
                    for (int ny = 0; ny < h4; ny++)
                    {
                        for (int nx = 0; nx < w4; nx++)
                        {
                            //ALEX FIX
                            //apply weight
                            int nPix4 = ny * w4 + nx;
                            int nPix2 = (ny / 2) * w2 + (nx / 2);
                            float wt = img1w[nPix4];
                            float v = img1[nPix4];
                            imgb[nPix2] += v;
                            imgbw[nPix2] += wt;
                        }
                    }

                    for (int ny = 0; ny < h2; ny++)
                    {
                        for (int nx = 0; nx < w2; nx++)
                        {
                            //ALEX FIX
                            //apply weight
                            int nPix2 = ny * w2 + nx;
                            float v = imgb[nPix2];
                            float t = imgbw[nPix2];
                            if (t == 0) { continue; }
                            //if (0 != t - Math.Floor(t))
                            //{
                            //    continue;
                            //}
                            t = 1 / t;
                            imgOut[nPix2] = v * t;
                        }
                    }
                    //ALEX
                    //Separable Filter
                    //ALEX icd 
                    //ALEX ! Frame Balance for Switching !
                    //TODO: write out answer
                    bw = new BinaryWriter(File.Create(@"slant.raw"));
                    for (int n = 0; n < imgOut.Length; n++)
                    {
                        bw.Write(imgOut[n]);
                    }
                    bw.Close();
                    bw.Dispose();
                    //  ezTiff.WriteTiff(@"C:\tmp\slant9.tif", imgOut, h2, w2);
                }
                //
                //Iterative "Improvement"
                //
                //ALEX Check for pixel center shift
                //
                // Back Project
                //
                //////var BP = new float[imgInp.pix.Length * h * y];
                ////////pixel box filter
                //////BoxSmooth(imgb, imgbw, img4PixBox, 4, 4);
                ////////motion blur
                //////InplaceBoxSmoothY(img4PixBox, Img4mBox, 0, 4);

                //////for (int nz = 0; nz < imgInp.pix.Length; nz++)
                //////{
                //////    for (int ny = 0; ny < h; ny++)
                //////    {
                //////        for (int nx = 0; nx < w; nx++)
                //////        {
                //////            //Add blank columns 
                //////            //- change the source index by dnx
                //////            int dnx = (nx / 129) + 1;
                //////            if (nx % 129 == 0)
                //////            {
                //////                continue;
                //////            }
                //////            //unbin
                //////            int nx4_0 = nx * 4;
                //////            int ny4_0 = ny * 4;
                //////            // Shear X for tilt
                //////            int nx4 = nx4_0 - ((ny4_0 * 80) / 1000);
                //////            //Motion in y per frame
                //////            int ny4 = ny4_0 - ((nz * 4 * 1019) / 1000);
                //////            ny4 += h4 - (4 * 256);

                //////            int niPix = ny * w + nx - dnx;
                //////            int nPix4 = ny4 * w4 + nx4;

                //////            float sample = img4mBox[nPix4];
                //////            BP[nz][niPix] = sample;
                //////            //calc data error
                //////            imgDiff[nz][niPix] =
                //////                imgInp.pix[nz][niPix] * imgGain[niPix]
                //////                - BP[nz][niPix];
                //////        }
                //////    }
                //////}
                //////for (int ny = 0; ny < h4; ny++)
                //////{
                //////    for (int nx = 0; nx < w4; nx++)
                //////    {
                //////        penalty[ny * w4 + nx] = localLaplacian(img1, img1w, ny, nx);
                //////    }
                //////}
                //
                //ALEX! DO implicit corrections between Tiles !
                // try svd of the 6000 frames
                //
                // reconstruct density?
                //
                // Forward project to model +weights
                //
                // pixel and motion blur
                // resample to Back project
                // calculate data fit diffence
                // calculate Laplacian Penalty
                //
                //Subtract Laplacian penalty
                // 
                //forward project sharpened diff
                //add
                // apply to img model
                //////////////////////////////////////////
                // refine fit of motion to data

                bool makeGain = false;
                if (makeGain)
                {
                    //float[,] svdData = new float[100, w * h];
                    //for (int nz = 0; nz < 100; nz++)
                    //{
                    //    int nframe = (nz * imgInp.pix.Length) / 100;
                    //    var imgn = imgInp.pix[nframe];
                    //    imgn.CopyTo(svdData, nz * w * h);
                    //}
                    //var dataFrame = new float[w * h];
                    //var svdDataSrc = OpenCvSharp.InputArray.Create(svdData, MatType.CV_32FC1);//[100, w * h];
                    //                                                                          //                var wd = OpenCvSharp.InputOutputArray.Create(new float[100, w * h], MatType.CV_32FC1);
                    //var u = new float[100, w * h];
                    //var vt = new float[100, w * h];
                    ////                svdDataSrc = new float[100, w * h];
                    ////                Cv2.SVDecomp(svdDataSrc, wd, u, vt, SVD.Flags.ModifyA);

                    var frameAvg0 = new float[w * h];
                    var nFrameAvg0 = new float[w * h];

                    for (int nz = 0; nz < imgpix.Length; nz++)
                    {
                        var imgn = imgpix[nz];
                        if (nz % 1000 == 0) { Console.Write("*"); }
                        for (int ny = 0; ny < h; ny++)
                        {
                            for (int nx = 0; nx < w; nx++)
                            {
                                int niPix = ny * w + nx;
                                frameAvg0[niPix] += imgn[niPix];
                                nFrameAvg0[niPix] += 1;
                            }
                        }
                    }
                    for (int ny = 0; ny < h; ny++)
                    {
                        for (int nx = 0; nx < w; nx++)
                        {
                            int nPix = ny * w + nx;
                            float t = Math.Max(1, nFrameAvg0[nPix]);
                            imgw[nPix] = 1.0f / t;
                            frameAvg0[nPix] /= t;
                        }
                    }

                    var frameAvg = new float[w * h];
                    var nFrameAvg = new float[w * h];
                    //Array.Clear(frameAvg, 0, frameAvg.Length);
                    //Array.Clear(nFrameAvg, 0, nFrameAvg.Length);
                    for (int nz = 0; nz < imgpix.Length; nz++)
                    {
                        float t = 0;
                        var imgn = imgpix[nz];
                        if (nz % 1000 == 0) { Console.Write("*"); }
                        for (int ny = 0; ny < h; ny++)
                        {
                            for (int nx = 0; nx < w; nx++)
                            {
                                int niPix = ny * w + nx;
                                var vPix = imgn[niPix];
                                //diff img to avg
                                t = (frameAvg0[niPix] - vPix) / (1 + frameAvg0[niPix]);
                                if (t > 0.5)
                                {
                                    //goto skip;
                                }
                                else
                                { t = t; }
                            }
                        }
                        for (int ny = 0; ny < h; ny++)
                        {
                            for (int nx = 0; nx < w; nx++)
                            {
                                int niPix = ny * w + nx;
                                var vPix = imgn[niPix];
                                frameAvg[niPix] += vPix;
                                nFrameAvg[niPix] += 1;
                            }
                        }
                    skip:;
                    }
                    for (int ny = 0; ny < h; ny++)
                    {
                        for (int nx = 0; nx < w; nx++)
                        {
                            int nPix = ny * w + nx;
                            float t = Math.Max(1, nFrameAvg[nPix]);
                            imgw[nPix] = 1.0f / t;
                            frameAvg[nPix] /= t;
                        }
                    }
                    //ezTiff.WriteTiff(fnGain + @".tiff", frameAvg, h, w);
                    //TODO: write out gain
                    //ezTiff.WriteTiff(@"C:\tmp\slantGainAvg.tif", frameAvg, h, w);
                }

                //Mat src = new Mat(@"c:\tmp\lenna.png", ImreadModes.Grayscale);
                // Mat src = Cv2.ImRead("lenna.png", ImreadModes.Grayscale);
                //Mat dst = new Mat();

                // Cv2.Canny(src, dst, 50, 200);
                //using (new Window("src image", src))
                //using (new Window("dst image", dst))
                //{
                //    Cv2.WaitKey();
                //}
                //Mat dst2 = new Mat(h4, w4, MatType.CV_8UC1);
                ////dst.SetArray(0,0,img0);
                ////dst = new Mat(40,256,MatType.CV_32FC1,&img0[0]);
                //using (new Window("dst image", dst2))
                //{
                //    Cv2.WaitKey();
                //}
            }
        }
    }
}
