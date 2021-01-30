import { async, ComponentFixture, TestBed } from '@angular/core/testing';

import { VcfSummaryComponent } from './vcf-summary.component';

describe('VcfSummaryComponent', () => {
  let component: VcfSummaryComponent;
  let fixture: ComponentFixture<VcfSummaryComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      declarations: [ VcfSummaryComponent ]
    })
    .compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(VcfSummaryComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
