#include <GL/gl.h>
#include <GL/glu.h>
#include <SDL2/SDL.h>
#include <iostream>

#define SDL_WINDOW_TITLE "SDL2"
#define SDL_WINDOW_WIDTH (640)
#define SDL_WINDOW_HEIGHT (480)

#define SDL_PrintError(name)                                   \
  do                                                           \
  {                                                            \
    std::cerr << #name << ": " << SDL_GetError() << std::endl; \
  } while (0)

static SDL_Window *gWindow;
static SDL_GLContext context;

static bool initialize()
{
  if (SDL_Init(SDL_INIT_EVERYTHING) < 0)
  {
    SDL_PrintError(SDL_Init);
    return false;
  }

  gWindow = SDL_CreateWindow(SDL_WINDOW_TITLE, SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
                             SDL_WINDOW_WIDTH, SDL_WINDOW_HEIGHT,
                             SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE);
  if (gWindow == nullptr)
  {
    SDL_PrintError(SDL_CreateWindow);
    goto err1;
  }

  // create OpenGL Context
  context = SDL_GL_CreateContext(gWindow);
  if (!context)
  {
    SDL_PrintError(SDL_GL_CreateContext);
    goto err2;
  }

  return true;

err2:
  SDL_DestroyWindow(gWindow);
err1:
  SDL_Quit();
  return false;
}

static bool initGL()
{

  /* Enable smooth shading */
  glShadeModel(GL_SMOOTH);

  /* Set the background black */
  glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

  /* Depth buffer setup */
  glClearDepth(1.0f);

  /* Enables Depth Testing */
  glEnable(GL_DEPTH_TEST);

  /* The Type Of Depth Test To Do */
  glDepthFunc(GL_LEQUAL);

  /* Really Nice Perspective Calculations */
  glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

  return (TRUE);
}

static bool resizeWindow(int width, int height)
{
  /* Height / width ration */
  GLfloat ratio;

  /* Protect against a divide by zero */
  if (height == 0)
    height = 1;

  ratio = (GLfloat)width / (GLfloat)height;

  /* Setup our viewport. */
  glViewport(0, 0, (GLsizei)width, (GLsizei)height);

  /* change to the projection matrix and set our viewing volume. */
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  /* Set our perspective */
  gluPerspective(45.0f, ratio, 0.1f, 100.0f);

  /* Make sure we're chaning the model view and not the projection */
  glMatrixMode(GL_MODELVIEW);

  /* Reset The View */
  glLoadIdentity();

  return (true);
}

static void finalize()
{
  SDL_GL_DeleteContext(context);
  SDL_DestroyWindow(gWindow);
  SDL_Quit();
}

static void render()
{
  /* Clear The Screen And The Depth Buffer */
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  /* Move Left 1.5 Units And Into The Screen 6.0 */
  glLoadIdentity();
  glTranslatef(-1.5f, 0.0f, -6.0f);

  glBegin(GL_TRIANGLES);          /* Drawing Using Triangles       */
  glColor3f(1.0f, 0.0f, 0.0f);    /* Red                           */
  glVertex3f(0.0f, 1.0f, 0.0f);   /* Top Of Triangle               */
  glColor3f(0.0f, 1.0f, 0.0f);    /* Green                         */
  glVertex3f(-1.0f, -1.0f, 0.0f); /* Left Of Triangle              */
  glColor3f(0.0f, 0.0f, 1.0f);    /* Blue                          */
  glVertex3f(1.0f, -1.0f, 0.0f);  /* Right Of Triangle             */
  glEnd();                        /* Finished Drawing The Triangle */

  /* Move Right 3 Units */
  glTranslatef(3.0f, 0.0f, 0.0f);

  /* Set The Color To Blue One Time Only */
  glColor3f(0.5f, 0.5f, 1.0f);

  glBegin(GL_QUADS);              /* Draw A Quad              */
  glVertex3f(1.0f, 1.0f, 0.0f);   /* Top Right Of The Quad    */
  glVertex3f(-1.0f, 1.0f, 0.0f);  /* Top Left Of The Quad     */
  glVertex3f(-1.0f, -1.0f, 0.0f); /* Bottom Left Of The Quad  */
  glVertex3f(1.0f, -1.0f, 0.0f);  /* Bottom Right Of The Quad */
  glEnd();                        /* Done Drawing The Quad    */

  /* Draw it to the screen */
  SDL_GL_SwapWindow(gWindow);

  /** NOTE: Sleep 10 msec. */
  SDL_Delay(10);
}

static bool input()
{
  SDL_Event event;

  /** NOTE: SDL_PollEvent does not sleep while SDL_WaitEvent sleep
      till event comes. SDL_WaitEvent is more relaxible than
      SDL_PollEvent. If input is combined with rendering,
      SDL_WaitEvent cannot be used. */
  while (SDL_PollEvent(&event))
  {
    switch (event.type)
    {
    case SDL_QUIT:
      return true;
      break;
    case SDL_WINDOWEVENT:
      switch (event.window.event)
      {
      case SDL_WINDOWEVENT_RESIZED:
        /* handle resize event */
        SDL_SetWindowSize(gWindow, event.window.data1, event.window.data2);
        resizeWindow(event.window.data1, event.window.data2);
        break;
      }
      break;
    default:
      break;
    }
  }

  return false;
}

int main(int argc, char *argv[])
{
  if (!initialize())
    return 1;

  initGL();
  resizeWindow(SDL_WINDOW_WIDTH, SDL_WINDOW_HEIGHT);

  while (1)
  {
    if (input())
      break;
    render();
  }

  finalize();
  return 0;
}