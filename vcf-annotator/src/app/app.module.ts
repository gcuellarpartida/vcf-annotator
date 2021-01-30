import { BrowserModule } from '@angular/platform-browser';
import { NgModule } from '@angular/core';

import { AppRoutingModule } from './app-routing.module';
import { AppComponent } from './app.component';
import { UploadBarComponent } from './upload-bar/upload-bar.component';
import { VcfSummaryComponent } from './vcf-summary/vcf-summary.component';
import { VcfAnnotationComponent } from './vcf-annotation/vcf-annotation.component';

@NgModule({
  declarations: [
    AppComponent,
    UploadBarComponent,
    VcfSummaryComponent,
    VcfAnnotationComponent
  ],
  imports: [
    BrowserModule,
    AppRoutingModule
  ],
  providers: [],
  bootstrap: [AppComponent]
})
export class AppModule { }
