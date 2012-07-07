//=============================================================================
//
//  CLASS GP_SmootherToolboxWidget
//
//  $Date$
//
//=============================================================================

#include "ui_SmootherToolbox.hh"
#include <QtGui>

class SmootherToolboxWidget : public QWidget, public Ui::PMC_SmootherToolbox
{
  Q_OBJECT

  public:
    SmootherToolboxWidget(QWidget *parent = 0);
};
