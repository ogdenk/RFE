#include "itkTestingExtractSliceImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryMedianImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkImageFileReader.h"
#include "itkScalarImageToTextureFeaturesFilter.h"
#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkBinaryImageToStatisticsLabelMapFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "QuickView.h"
#include <stdio.h>
#include <vector>


// These are the Freesurfer label values for the structures of interest
#define LeftHippocampus 17
#define RightHippocampus 53
#define LeftThalamus 10
#define RightThalamus 49


int main(int argc, char * argv[])
{
  
  unsigned short int labelValue; // holds the current label value during the loop
  std::string labelName;		 // holds the current label name during the loop
  const unsigned int Dimension = 3;
  int i, j, k; // used for general purpose indexing
  int radius=1;  // size of neighborhood for GLCM determination
  bool found;  // used to build the texture neighborhood vector

  // both image data and labels can be opened as UINT's
  typedef  unsigned short  int								   PixelType;
  typedef itk::Image< PixelType, Dimension >                   ImageType;
  typedef itk::Image< PixelType, Dimension-1>				   SliceType;

  ImageType::SizeType medianRadius; // used for binary median image filter
  medianRadius[0] = 1;
  medianRadius[1] = 1;
  medianRadius[2] = 1;

  // Get the input directory specified on the command line
  std::string inputDirectory;
  if (argc > 1)
  {
	 	inputDirectory = argv[1];
  }
  else
  {
  	std::cout <<  "Input a directory containing files to process" << std::endl;
	exit(1);
  }

  // the file names for the brain scan and the label image
  std::string imageData = inputDirectory + "/orig.nrrd";
  std::string labelData = inputDirectory + "/labels.nrrd";
  std::string outputstring; // used to format text output to the file

  // Open the output file and put the directory name on the first line
  FILE * outputFile = fopen((inputDirectory + "features.csv").c_str(), "wt");
  outputstring = "Output Directory = " + inputDirectory + "\n";
  fprintf(outputFile, outputstring.c_str());
	
  // Create vectors of label names and label values
  std::vector<std::string> labelNames;
  labelNames.reserve(4); //storage for four labels
  labelNames.push_back("LeftHippocampus");
  labelNames.push_back("RightHippocampus");
  labelNames.push_back("LeftThalamus");
  labelNames.push_back("RightThalamus");

  std::vector<unsigned short int> labelValues;
  labelValues.reserve(4);  //storage for four values
  labelValues.push_back(LeftHippocampus);
  labelValues.push_back(RightHippocampus);
  labelValues.push_back(LeftThalamus);
  labelValues.push_back(RightThalamus);

  // Create iterators for the vectors.  We will iterate over the structures when doing the processing
  std::vector<std::string>::iterator labelNamesIterator = labelNames.begin();
  std::vector<unsigned short int>::iterator labelValuesIterator = labelValues.begin();

  // open the image data
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  ReaderType::Pointer imageReader = ReaderType::New();
  imageReader->SetFileName(imageData );
  imageReader->Update();

  // open the labels
  typedef itk::ImageFileReader< ImageType >  LabelReaderType;
  LabelReaderType::Pointer labelReader = LabelReaderType::New();
  labelReader->SetFileName(labelData);
  labelReader->Update();
  labelReader->GetOutput()->SetSpacing(imageReader->GetOutput()->GetSpacing());
  labelReader->GetOutput()->SetOrigin(imageReader->GetOutput()->GetOrigin());

 /* // This is just to get a look at one of the slices of the image data.  This block can be commented out
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ImageType::SizeType imagesize = labelReader->GetOutput()->GetLargestPossibleRegion().GetSize();
  ImageType::RegionType extractRegion;
  extractRegion.SetSize(0,imagesize[0]);
  extractRegion.SetSize(1,imagesize[1]);
  extractRegion.SetSize(2,0);
  extractRegion.SetIndex(0,0);
  extractRegion.SetIndex(1,0);
  extractRegion.SetIndex(2,128); // extract image slice 128 for viewing purposes

  
  typedef itk::Testing::ExtractSliceImageFilter<ImageType, SliceType> ExtractSliceType;
  ExtractSliceType::Pointer SliceExtracter = ExtractSliceType::New();
  SliceExtracter->SetInput(labelReader->GetOutput());
  SliceExtracter->SetExtractionRegion(extractRegion);
  SliceExtracter->SetDirectionCollapseToIdentity();
  SliceExtracter->Update();

  // view a slice
  QuickView viewer;
  viewer.AddImage(SliceExtracter->GetOutput());
  viewer.Visualize();
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  */

  // Create smoothed version of the image data
  typedef itk::DiscreteGaussianImageFilter<ImageType, ImageType> SmoothingType;
  SmoothingType::Pointer smoother = SmoothingType::New();
  smoother->SetInput(imageReader->GetOutput());
  smoother->SetVariance(1.0);
  smoother->SetUseImageSpacingOn(); // this is the default value but just in case
  smoother->SetMaximumKernelWidth(5);

  // Set up the pipeline then we will iterate over the labels to be processed
  /////////////////////////////////////////////////////////////////////////////////////////

  //  Median filter to smooth the edges of the labels
  typedef itk::BinaryMedianImageFilter<ImageType, ImageType> MedianFilterType;
  MedianFilterType::Pointer medianFilter = MedianFilterType::New();
  medianFilter->SetRadius(medianRadius);
  medianFilter->SetBackgroundValue(0);
  medianFilter->SetInput(labelReader->GetOutput());

  //  Shape filter to extract shape metrics from the labels
  typedef itk::BinaryImageToShapeLabelMapFilter<ImageType> BinaryImageToShapeLabelMapFilterType;
  BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToShapeLabelMapFilter = BinaryImageToShapeLabelMapFilterType::New();
  binaryImageToShapeLabelMapFilter->SetComputeFeretDiameter(true);
  binaryImageToShapeLabelMapFilter->SetInput(medianFilter->GetOutput());
  BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType* labelObject;

  //  Statistics filter to get the gray level statistics from the MRI image for a particular label location
  typedef itk::BinaryImageToStatisticsLabelMapFilter<ImageType, ImageType> BinaryImageToStatisticsLabelMapFilterType;
  BinaryImageToStatisticsLabelMapFilterType::Pointer BinaryToStatisticsFilter = BinaryImageToStatisticsLabelMapFilterType::New();
  BinaryToStatisticsFilter->SetInput1(medianFilter->GetOutput());
  BinaryToStatisticsFilter->SetComputeHistogram(TRUE);
  BinaryToStatisticsFilter->SetCoordinateTolerance(0.01);
  BinaryImageToStatisticsLabelMapFilterType::OutputImageType::LabelObjectType* StatlabelObject;
 
  //  Texture filter to get the Haralick features
  typedef itk::Statistics::ScalarImageToTextureFeaturesFilter<ImageType> TextureFilterType;
  TextureFilterType::Pointer textureFilter = TextureFilterType::New();
  textureFilter->SetMaskImage(medianFilter->GetOutput());
  textureFilter->SetFastCalculations(false);
  textureFilter->SetNumberOfBinsPerAxis(1024);
  textureFilter->SetPixelValueMinMax(0, 2048);


  TextureFilterType::OffsetType   offset, test;
  TextureFilterType::OffsetVectorPointer   offset1;
  offset1 = TextureFilterType::OffsetVector::New();
  TextureFilterType::OffsetVector::ConstIterator vIt; 
  const TextureFilterType::FeatureValueVector* output;
  const TextureFilterType::FeatureValueVector* outputSD;

  // Build the offset vector for the texture filter
  radius = 1;
  for (i = -radius; i <= radius; i++)
  {
	  for (j = -radius; j <= radius; j++)
	  {
		  for (k = -radius; k <= radius; k++)
		  {
			  if (((i + j + k) <= 0) && !(i == 0 && j == 0 && k == 0))
			  {
				  found = false;
				  offset[0] = i;
				  offset[1] = j;
				  offset[2] = k;
				  test[0] = -i;
				  test[1] = -j;
				  test[2] = -k;
				  for (vIt = offset1->Begin(); vIt != offset1->End(); ++vIt)
				  {
					  if (vIt.Value() == test)
						  found = true;
				  }

				  if (!found)
					  offset1->push_back(offset);
			  }
		  }
	  }
  }
  // Now set the offset vector for the texture filter
  textureFilter->SetOffsets(offset1);


  // Iterate over the labels we want to process
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  i = 0;
  while (labelNamesIterator != labelNames.end())
  {
	  labelName = labelNames[i];
	  labelValue = labelValues[i];

	  // First we want to extract the current label and median filter to remove rough edges.
	  // need to set the appropriate value in the median filter
	  medianFilter->SetForegroundValue(labelValue);
	  medianFilter->Update();

	  // Now get the shape statistics for that label
	  binaryImageToShapeLabelMapFilter->SetInputForegroundValue(labelValue);
	  binaryImageToShapeLabelMapFilter->Update();
	  std::cout << "There is " << binaryImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() << " object." << std::endl;

	  // Loop over the regions (should only be 1)
	  std::cout << "Shape values for " << labelName << std::endl;
	  std::cout << "Roundness, Flatness, Feret Diameter, Number of Pixels, Physical Size " << std::endl;
	  for (j = 0; j < binaryImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects(); j++)
	  {
		  labelObject = binaryImageToShapeLabelMapFilter->GetOutput()->GetNthLabelObject(j);
		  std::cout << labelObject->GetRoundness() << std::endl;
		  outputstring = labelName + "Roundness, %f \n";
		  fprintf(outputFile, outputstring.c_str(), labelObject->GetRoundness());
		  std::cout << labelObject->GetFlatness() << std::endl;
		  outputstring = labelName + "Flatness, %f \n";
		  fprintf(outputFile, outputstring.c_str(), labelObject->GetFlatness());
		  std::cout << labelObject->GetFeretDiameter() << std::endl;
		  outputstring = labelName + "FeretDiameter, %f \n";
		  fprintf(outputFile, outputstring.c_str(), labelObject->GetFeretDiameter());
		  std::cout << labelObject->GetPhysicalSize() << std::endl;
		  outputstring = labelName + "Volume, %f \n";
		  fprintf(outputFile, outputstring.c_str(), labelObject->GetPhysicalSize());
		  std::cout << labelObject->GetElongation() << std::endl;
		  outputstring = labelName + "Elongation, %f \n";
		  fprintf(outputFile, outputstring.c_str(), labelObject->GetElongation());

		  //std::cout << "Object " << i << " has principal axes " << labelObject->GetPrincipalAxes() << std::endl;
		  //std::cout << "Object " << i << " has principal moments " << labelObject->GetPrincipalMoments() << std::endl;
		  
		  std::cout << std::endl << std::endl;
	  }

	  // Now get the statistics for the un-Blurred image data using the current label

	  BinaryToStatisticsFilter->SetInput2(imageReader->GetOutput());
	  BinaryToStatisticsFilter->SetInputForegroundValue(labelValue);
	  BinaryToStatisticsFilter->Update();

	  std::cout << "There is " << BinaryToStatisticsFilter->GetOutput()->GetNumberOfLabelObjects() << " object with statistics." << std::endl;
	  std::cout << "Statistics values from un-blurred image for " << labelName << std::endl;
	  std::cout << "Mean, Median, Skewness, Kurtosis, Standard Deviation " << std::endl;
	  for (k = 0; k < BinaryToStatisticsFilter->GetOutput()->GetNumberOfLabelObjects(); k++)
	  {
		  StatlabelObject = BinaryToStatisticsFilter->GetOutput()->GetNthLabelObject(k);
		  // Output the shape properties of the ith region
		  // Feret diameter, mean, median, skewness, kurtosis, sigma
		  std::cout << StatlabelObject->GetMean() << std::endl;
		  std::cout << StatlabelObject->GetMedian() << std::endl;
		  std::cout << StatlabelObject->GetSkewness() << std::endl;
		  std::cout << StatlabelObject->GetKurtosis() << std::endl;
		  std::cout << StatlabelObject->GetStandardDeviation() << std::endl;
		  std::cout << std::endl << std::endl;
	  }

	  // Now get the texture features for the unblurred image
	  textureFilter->SetInput(imageReader->GetOutput());
	  textureFilter->SetInsidePixelValue(labelValue);
	  textureFilter->Update();

	  output = textureFilter->GetFeatureMeans();
	  outputSD = textureFilter->GetFeatureStandardDeviations();

	  std::cout << "Radius 1 Texture Features for: " << labelName << std::endl;
	  for (unsigned int i = 0; i < output->size(); ++i)
	  {
		  std::cout << (*output)[i] << std::endl;
		  std::cout << (*outputSD)[i] << std::endl;
	  }
	  std::cout << std::endl << std::endl;

	  // Now get the statistics for the Blurred image data using the current label

	  BinaryToStatisticsFilter->SetInput2(smoother->GetOutput());
	  BinaryToStatisticsFilter->SetInputForegroundValue(labelValue);
	  BinaryToStatisticsFilter->Update();

	  std::cout << "There is " << BinaryToStatisticsFilter->GetOutput()->GetNumberOfLabelObjects() << " object with statistics." << std::endl;
	  std::cout << "Statistics values from blurred image for " << labelName << std::endl;
	  std::cout << "Mean, Median, Skewness, Kurtosis, Standard Deviation " << std::endl;
	  for (k = 0; k < BinaryToStatisticsFilter->GetOutput()->GetNumberOfLabelObjects(); k++)
	  {
		  StatlabelObject = BinaryToStatisticsFilter->GetOutput()->GetNthLabelObject(k);
		  // Output the shape properties of the ith region
		  // Feret diameter, mean, median, skewness, kurtosis, sigma
		  std::cout << StatlabelObject->GetMean() << std::endl;
		  std::cout << StatlabelObject->GetMedian() << std::endl;
		  std::cout << StatlabelObject->GetSkewness() << std::endl;
		  std::cout << StatlabelObject->GetKurtosis() << std::endl;
		  std::cout << StatlabelObject->GetStandardDeviation() << std::endl;
		  std::cout << std::endl << std::endl;
	  }

	  // Now get the texture features for the blurred image
	  textureFilter->SetInput(smoother->GetOutput());
	  textureFilter->SetInsidePixelValue(labelValue);
	  textureFilter->Update();

	  output = textureFilter->GetFeatureMeans();
	  outputSD = textureFilter->GetFeatureStandardDeviations();

	  std::cout << "Radius 1 Texture Features (blurred image) for: " << labelName << std::endl;
	  for (unsigned int i = 0; i < output->size(); ++i)
	  {
		  std::cout << (*output)[i] << std::endl;
		  std::cout << (*outputSD)[i] << std::endl;
	  }




	  i++;
	  labelValuesIterator++;
	  labelNamesIterator++;
  }

  fclose(outputFile);
 
  std::cout << std::endl;
  std::cout << std::endl;


  getchar();
}