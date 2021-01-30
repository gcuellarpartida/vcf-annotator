import { async, ComponentFixture, TestBed } from '@angular/core/testing';

import { VcfAnnotationComponent } from './vcf-annotation.component';

describe('VcfAnnotationComponent', () => {
  let component: VcfAnnotationComponent;
  let fixture: ComponentFixture<VcfAnnotationComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      declarations: [ VcfAnnotationComponent ]
    })
    .compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(VcfAnnotationComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
