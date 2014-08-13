#include <vtkVersion.h>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleImage.h>
#include <vtkRenderer.h>
#include <vtkImageMapper.h>
#include <vtkActor2D.h>
#include <vtkMath.h>

static void CreateColorImage(vtkImageData*);

int main(int, char *[])
{
  vtkSmartPointer<vtkImageData> colorImage = vtkSmartPointer<vtkImageData>::New();
  CreateColorImage(colorImage);
  
  vtkSmartPointer<vtkImageMapper> imageMapper = vtkSmartPointer<vtkImageMapper>::New();
#if VTK_MAJOR_VERSION <= 5
  imageMapper->SetInputConnection(colorImage->GetProducerPort());
#else
  imageMapper->SetInputData(colorImage);
#endif
  imageMapper->SetColorWindow(255);
  imageMapper->SetColorLevel(127.5);
  
  vtkSmartPointer<vtkActor2D> imageActor = vtkSmartPointer<vtkActor2D>::New();
  imageActor->SetMapper(imageMapper);
  imageActor->SetPosition(400, 400);
        
  // Setup renderers
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  
  // Setup render window
  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  
  renderWindow->AddRenderer(renderer);

  // Setup render window interactor
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    
  vtkSmartPointer<vtkInteractorStyleImage> style = vtkSmartPointer<vtkInteractorStyleImage>::New();
  
 // renderWindowInteractor->SetInteractorStyle(style);

  // Render and start interaction
  renderWindowInteractor->SetRenderWindow(renderWindow);
  
  //renderer->AddViewProp(imageActor);
  renderer->AddActor2D(imageActor);

  renderWindow->Render();
  renderer->ResetCamera();
  renderWindow->SetSize(1200,1200); //(width, height)
  renderWindowInteractor->Start();
  
  return EXIT_SUCCESS;
}

void CreateColorImage(vtkImageData* image)
{
  unsigned int dim = 20;
  unsigned int block = 20;
  
  image->SetDimensions(dim*block, dim*block, 1);
#if VTK_MAJOR_VERSION <= 5
  image->SetNumberOfScalarComponents(3);
  image->SetScalarTypeToUnsignedChar();
  image->AllocateScalars();
#else
  image->AllocateScalars(VTK_UNSIGNED_CHAR,3);
#endif
  for(unsigned int x = 0; x < dim; x++)
      for(unsigned int y = 0; y < dim; y++)
      {
          double pixelval= vtkMath::Random(0.0, 255.0);
          for(unsigned int bx = 0; bx < block; bx++)
              for(unsigned int by = 0; by < block; by++)
              {
                  double* pixel = static_cast<double*>(image->GetScalarPointer(x*block+bx,y*block+by,0));
                  if (x==0)
                      pixel[0] = pixelval;
              }
      }
   
    
  image->Modified();
}
